/*
 * controller_jair.c
 *
 * Implementacion fiel al articulo COMRob 2025:
 *   - PID posicion  (ec.12 + 11a)
 *   - STSMC orientacion (ec.18 + 11b)
 *   - Observador SM-FTO (ec.5) para xi y eta
 *
 * ================================================================
 *  NOTAS CRITICAS DE IMPLEMENTACION
 * ================================================================
 *
 * 1. ESCALA DE TORQUES (ec.11b):
 *    El paper escribe tau = J^-1*(tau_bar - ...).
 *    El firmware Crazyflie espera control->torque[] en [N*m].
 *    power_distribution() aplica J^-1 internamente.
 *    => Se entrega tau_bar directamente en [N*m] SIN multiplicar J^-1.
 *    => Las ganancias k1, k2 del paper (escala rad/s^2) se escalan por J:
 *         k1_Nm = k1_paper * J  ~ 10.1181 * 9.827e-5 ~ 9.94e-4
 *         k2_Nm = k2_paper * J  ~ 50.05   * 9.827e-5 ~ 4.92e-3
 *    => k0 NO se escala (es adimensional, superficie en rad/s).
 *
 * 2. ERROR DE YAW (causa del giro continuo):
 *    El error angular debe estar en [-pi, pi].
 *    wrap_angle() normaliza cualquier diferencia angular.
 *
 * 3. OBSERVADOR (ec.5):
 *    Se integra con Euler hacia adelante a la misma tasa que
 *    el controlador (ATTITUDE_RATE).
 *    El observador converge en tiempo finito T0; antes de T0 las
 *    estimaciones son ruidosas. Se expone obs_enabled para activar
 *    el observador solo cuando el dron este estable.
 *    xi_hat3 y eta_hat3 se usan en ec.11a y 11b respectivamente.
 *
 * 4. FUNCIONES DEL OBSERVADOR (distintas a las del STSMC):
 *    phi1(s) = ceil(s)^(2/3)  para estados .1 y .4
 *    phi2(s) = ceil(s)^(1/3)  para estados .2 y .5 (NO 1/2)
 *    phi3(s) = sign(s)        para estados .3 y .6
 *    Estas son las funciones exactas del paper (Seccion IV).
 *
 * 5. CAPA LIMITE:
 *    Para el STSMC se usan capas limite (delta_s, delta_e) para
 *    eliminar chattering. Para el observador se usa sign() puro
 *    porque phi3 es sign() y el observador necesita la ganancia
 *    maxima para la identificacion de perturbaciones.
 */

#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "math3d.h"
#include "controller_jair.h"
#include "physicalConstants.h"
#include "power_distribution.h"
#include "platform_defaults.h"

/* ============================================================== */
/*  Constantes físicas del CF2.1 (platform_defaults_cf2.h)       */
/* ============================================================== */
/*
 * Derivadas de platform_defaults_cf2.h:
 *   THRUST_MAX    = 0.12 N/motor
 *   THRUST_MIN    = 0.012818 N/motor
 *   ARM_LENGTH    = 0.046 m
 *   THRUST2TORQUE = 0.006993 m
 *   CF_MASS       = 0.029 kg
 *
 * Torque máximo físico por eje:
 *   roll/pitch: (THRUST_MAX - THRUST_MIN) * 2 * ARM = 0.00986 N·m
 *   yaw:        (THRUST_MAX - THRUST_MIN) * 2 * T2T = 0.00150 N·m
 *
 * Se usa el 80% como límite de saturación para dejar margen.
 */
#define CF2_MAX_THRUST_N     0.480f   /* 4 * THRUST_MAX [N]           */
#define MAX_TORQUE_ROLL_NM   0.0079f  /* 80% de 0.00986 N·m           */
#define MAX_TORQUE_YAW_NM    0.0012f  /* 80% de 0.00150 N·m           */

/* ============================================================== */
/*  Ganancias por defecto (Tablas IV-VI del paper)                */
/* ============================================================== */
/*
 * Ganancias observador: valores exactos del paper (Seccion V)
 * K1..K3: subsistema lineal (posicion)
 * K4..K6: subsistema angular (orientacion)
 *
 * Ganancias STSMC phi/theta (Tabla IV-VI, fila STSMC):
 *   k0_phi = k0_theta = 2.6  [adimensional, superficie en rad/s]
 *   k0_psi = 2.1             [yaw usa ganancia diferente]
 *   k1_phi = k1_theta = 10.1181 * J_phi = 10.1181 * 9.827e-5 = 9.94e-4 N*m
 *   k1_psi = 4.3732 * J_psi = 4.3732 * 9.613e-5 = 4.20e-4 N*m
 *   k2_phi = k2_theta = 50.05 * J_phi = 50.05 * 9.827e-5 = 4.92e-3 N*m/s
 *   k2_psi = 9.35 * J_psi = 9.35 * 9.613e-5 = 8.99e-4 N*m/s
 */
static controllerJair_t g_self = {
  .mass = 0.032f,   /* Masa medida experimentalmente (hover 20Mar26) *
                   * Error Z=11.6cm implica masa_real~31.8g.       *
                   * Ajustar ctrlJair.mass desde cfclient si drift  *
                   * en Z persiste. Subir si sigue bajo, bajar      *
                   * si sube de mas.                               */
  //.J    = {9.827e-5f, 8.185e-5f, 9.613e-5f}, // Para CF 2.1
  .J = {16.571710e-6f, 16.655602e-6f, 29.261652e-6f}, //Para CF 2.0 (Lee)

  /* Posicion PID — ganancias para EXPERIMENTO REAL
   *
   * Las ganancias del paper (k1xi=28.5) son para simulacion.
   * En experimento real, un error de 30cm con k1xi=28.5 genera
   * phi*=41 deg -> volteo inmediato del quad.
   *
   * Regla: nu_x_max = Kpos_P * Kpos_P_limit < g*sin(20deg) = 3.35 m/s2
   * Con Kpos_P=2.0 y limit=0.5m: nu_x_max = 1.0 m/s2 -> phi*=5.9 deg
   *
   * Ajustar Kpos_P hacia arriba gradualmente en vuelo si la
   * respuesta es demasiado lenta.                                */
  /* Ganancias de posicion para EXPERIMENTO REAL:
   * Lee usa Kpos_P=7.0 con controlModeForceTorque exitosamente.
   * max_tilt=20deg garantiza no-volteo con cualquier error.
   *
   * Kpos_P=5.0: e_x=0.3m -> nu_x=1.5 m/s² -> phi*=8.8° (seguro)
   * Kpos_P=7.0: e_x=0.3m -> nu_x=2.1 m/s² -> phi*=12.4° (seguro)
   *
   * Kpos_D: amortiguamiento. Lee usa 4.0. Mantener igual.
   * Kpos_I: integral pequeña para corregir deriva lenta.
   *         Con I pequeño la deriva XY que se observó se corrige. */
  /* Ganancias ajustadas tras experimento de hover (20 Mar 2026):
   * - Error Z estacionario = 11.6cm → masa real ~31.8g
   *   Kpos_I_limit subido de 0.5 a 2.0 para compensar
   *   Kpos_I_z subido a 1.5 para convergencia mas rapida en Z
   * - Deriva XY inicial: Kpos_P_xy subido de 6 a 8 para
   *   recuperar mas rapido la posicion XY
   * - Masa actualizada a 0.032 (medida implicita del vuelo)
   *
   * Regla de seguridad (max_tilt=20deg):
   *   nu_x_max = Kpos_P_xy * Kpos_P_limit
   *   Con Kpos_P=8, limit=0.5: nu_x_max=4.0 m/s2 -> phi*=22deg (limite)
   *   => reducir Kpos_P_limit a 0.4 para mantener phi*<20deg          */
  .Kpos_P       = {8.0f,  8.0f,  6.0f},  /* xy mas agresivo, z igual  */
  .Kpos_D       = {4.0f,  4.0f,  4.0f},
  .Kpos_I       = {0.5f,  0.5f,  1.5f},  /* z con mas integral        */
  .Kpos_P_limit = 0.4f,   /* max 0.4m: phi*_max=18.8deg con Kpos_P=8  */
  .Kpos_D_limit = 1.0f,
  .Kpos_I_limit = 2.0f,   /* subido de 0.5: integral puede compensar Z */
  .max_tilt_rad = 0.3491f,  /* 20 deg — limite anti-volteo  */
  .Lambda_xi    = {0.0f, 0.0f, 0.0f},
  //.Lambda_xi    = {2.47e-7f, 2.47e-7f, 2.47e-7f},

  /* STSMC orientacion: ganancias escaladas a [N*m]
   *
   * Problema raíz del chattering:
   *   Ruido giroscopio CF2.1 ~ 0.05 rad/s en eps_eta.
   *   Con k0=2.6: contribución de k0*e^(2/3) con e=5deg=0.087rad
   *     => 2.6 * 0.087^(2/3) = 2.6 * 0.198 = 0.515 rad/s en s
   *   Eso hace que |s| >> delta_s=0.15 siempre, el controlador
   *   opera fuera de la capa límite y conmuta continuamente.
   *
   * Solución: reducir k0 para que con errores normales de vuelo
   *   (e~5deg, eps~0.1 rad/s), |s| < delta_s:
   *   s = eps_eta + k0*e^(2/3) ~ 0.1 + k0*0.2
   *   Para |s| < 0.15: k0 < (0.15-0.1)/0.2 = 0.25
   *   => usar k0 = 0.1 (muy conservador para primer vuelo)
   *
   * k1 escalado: k1_Nm = k1_paper * J
   *   phi/theta: 10.1181 * 9.827e-5 = 9.94e-4 N*m
   *   psi: reducido a 1/5 del papel para yaw más suave
   *
   * k2: integrador lento para evitar wind-up
   *   max_v por eje: psi tiene 5x menos que phi/theta
   */
  .Ksmc_k0    = {0.10f,   0.10f,   0.90f},
  .Ksmc_k1    = {9.94e-4f, 9.94e-4f, 4.20e-4f},
  //.Ksmc_k1    = {1.677e-4f, 1.677e-4f, 1.280e-4f},
  .Ksmc_k2    = {4.92e-3f, 4.92e-3f, 8.99e-4f},
  //.Ksmc_k2    = {8.294e-4f, 8.294e-4f, 2.736e-4f},
  .Ksmc_max_v  = {1.479e-3f, 1.479e-3f, 2.25e-4f},  /* 15% tau_max fisico */
  /* delta_s > ruido_gyro + k0*e_tipico para operar dentro de capa */
  .Ksmc_delta_s = 0.20f,   /* [rad/s]: > 0.05(ruido) + 0.1*0.2(k0*e) */
  .Ksmc_delta_e = 0.10f,   /* [rad]: ~5.7deg                          */
  .Lambda_eta = {0.0f, 0.0f, 0.0f},
  //.Lambda_eta = {1.52e-7f, 1.52e-7f, 1.52e-7f},
  .Xi_coef    = {0.0f, 0.0f, 0.0f}, 
  //.Xi_coef    = {-0.1453f, -0.0261f, 0.1708f}, // Para CF 2.1
  //.Xi_coef    = {-0.7604f, 0.7617f, -0.0031f}, //Para CF 2.0 (Lee)
  
  /* Observador SM-FTO — ganancias calculadas desde D
   *
   * CRITERIO (Teorema 1 del paper):
   *   K3 = D_xi  * obs_margin  (K3 > D_xi  garantiza convergencia)
   *   K6 = D_eta * obs_margin  (K6 > D_eta garantiza convergencia)
   * K1,K2,K4,K5: velocidad de convergencia (no dependen de D).
   *
   * El usuario configura D_xi y D_eta en cfclient según la
   * perturbación esperada. K3 y K6 se recalculan automáticamente.
   *
   * Tabla de referencia (turbulencia 30% a 3 rad/s):
   *   Sin perturb (lab):  D_xi=0.005, D_eta=0.0001
   *   Ventilador 0.84m/s: D_xi=0.061, D_eta=0.000039
   *   Soplido ~1.5m/s:    D_xi=0.194, D_eta=0.00013
   *   Soplido fuerte 3m/s: D_xi=0.775, D_eta=0.00052
   *
   * ACTIVACIÓN: obs_enabled=1 desde el inicio del vuelo es seguro
   * con D_xi conservador. El observador arranca con xi_hat=medición
   * (error cero) y converge suavemente.                           */
  /* obs_D_xi=0 → K3=0: observador solo estima xi_hat1 y xi_hat2
   * (sin acumular xi_hat3). Cambiar a 0.005+ para estimar perturbaciones. */
  //.obs_D_xi     = 0.0f,    /* 0 = solo velocidad/posición, sin perturbación */
  //.obs_D_eta    = 0.0f,    /* 0 = solo vel angular, sin perturbación   */
  //.obs_margin   = 1.20f,   /* K = D * 1.20 (20% sobre la cota)         */
  /* Observador SM-FTO (Seccion V del paper) */
  .Kobs_1 = {1.0766f, 1.0766f, 1.2324f},
  .Kobs_2 = {0.5925f, 0.5925f, 0.7256f},
  /* K3: ESCENARIO A (sin perturb). Cambiar a 0.0729 con ventilador */
  .Kobs_3 = {0.0060f, 0.0060f, 0.0060f},
  //.Kobs_3 = {0.1716f, 0.1716f, 0.2574f},
  .Kobs_4 = {1.4930f, 1.6083f, 1.6602f},
  .Kobs_5 = {0.9675f, 1.0817f, 1.1345f},
  /* K6: ESCENARIO A (sin perturb). Cambiar a 0.000047 con ventilador */
  .Kobs_6 = {0.0001f, 0.0001f, 0.0001f},
  //.Kobs_6 = {0.4576f, 0.5720f, 0.6292f},
  .obs_enabled   = 0,      /* activo desde el inicio del vuelo         */
  .obs_delta_xi  = 0.05f,  /* capa límite: filtrar ruido sin perturb   */
  .obs_delta_eta = 0.05f,  /* reducir a 0.05 cuando haya perturb real  */

  /* Filtro IIR para eta_dot_d: fc=5Hz @ 1kHz => alpha=0.969     */
  .eta2d_filter_alpha = 0.969f,
};

/* ============================================================== */
/*  Funciones auxiliares                                          */
/* ============================================================== */

static inline float clampf(float v, float lo, float hi)
{
  return v < lo ? lo : (v > hi ? hi : v);
}

static inline struct vec vclampscl(struct vec v, float lo, float hi)
{
  return mkvec(clampf(v.x, lo, hi),
               clampf(v.y, lo, hi),
               clampf(v.z, lo, hi));
}

/*
 * wrap_angle(a): normaliza un angulo a [-pi, pi].
 * Imprescindible para el error de yaw.
 */
static inline float wrap_angle(float a)
{
  while (a >  (float)M_PI) a -= 2.0f * (float)M_PI;
  while (a < -(float)M_PI) a += 2.0f * (float)M_PI;
  return a;
}

/*
 * ceil_pow(a, gamma): ceil(a)^gamma = |a|^gamma * sign(a)
 * Funcion basica de SMC, sin suavizado.
 * Se usa en el OBSERVADOR donde se necesita sign() puro (phi3).
 */
static inline float ceil_pow(float a, float gamma)
{
  if (fabsf(a) < 1e-12f) return 0.0f;
  if (gamma == 0.0f) return (a > 0.0f ? 1.0f : -1.0f);
  return powf(fabsf(a), gamma) * (a > 0.0f ? 1.0f : -1.0f);
}

static inline struct vec vceil_pow(struct vec v, float gamma)
{
  return mkvec(ceil_pow(v.x, gamma),
               ceil_pow(v.y, gamma),
               ceil_pow(v.z, gamma));
}

/*
 * ceil_pow_smooth(a, gamma, delta): igual que ceil_pow pero con
 * capa limite de ancho delta para eliminar chattering.
 * Se usa en el STSMC (no en el observador).
 *
 *   |a| >= delta : |a|^gamma * sign(a)
 *   |a| <  delta : |a|^gamma * (a/delta)   [continuo]
 *
 * Para gamma=0: equivale a sat(a/delta).
 */
static inline float ceil_pow_smooth(float a, float gamma, float delta)
{
  if (fabsf(a) < 1e-12f) return 0.0f;
  float base = (gamma == 0.0f) ? 1.0f : powf(fabsf(a), gamma);
  if (fabsf(a) >= delta) {
    return base * (a > 0.0f ? 1.0f : -1.0f);
  }
  return base * (a / delta);
}

static inline struct vec vceil_pow_smooth(struct vec v,
                                           float gamma, float delta)
{
  return mkvec(ceil_pow_smooth(v.x, gamma, delta),
               ceil_pow_smooth(v.y, gamma, delta),
               ceil_pow_smooth(v.z, gamma, delta));
}

/*
 * Termino de Coriolis: Xi * wη(η2)
 * wη(η2) = (θ̇ψ̇, φ̇ψ̇, φ̇θ̇)^T
 * Xi_coef = (bphi, btheta, bpsi) [adimensional]
 * Resultado en las mismas unidades que tau_bar [N*m] si
 * Xi_coef esta escalado por J.
 */
static inline struct vec coriolis_term(struct vec eta2, struct vec xi)
{
  return mkvec(xi.x * eta2.y * eta2.z,
               xi.y * eta2.x * eta2.z,
               xi.z * eta2.x * eta2.y);
}

/* Macro para copiar struct vec a campos log float escalares */
#define LOG_VEC3(self, prefix, vec) \
  do { \
    (self)->log_##prefix##_x = (vec).x; \
    (self)->log_##prefix##_y = (vec).y; \
    (self)->log_##prefix##_z = (vec).z; \
  } while(0)

/* ============================================================== */
/*  API publica                                                   */
/* ============================================================== */

void controllerJairReset(controllerJair_t *self)
{
  /* Controlador */
  self->i_error_pos  = vzero();
  self->v_smc        = vzero();
  self->prev_eta1_d  = vzero();
  self->eta2d_filt   = vzero();
  self->prev_xi_d2   = vzero();
  self->prev_valid      = false;

  /* Observador: inicializar con medicion actual si estuviera
   * disponible. Al reset, se pone a cero; convergera rapidamente. */
  self->xi_hat1  = vzero();
  self->xi_hat2  = vzero();
  self->xi_hat3  = vzero();
  self->eta_hat1 = vzero();
  self->eta_hat2 = vzero();
  self->eta_hat3 = vzero();
}

void controllerJairInit(controllerJair_t *self)
{
  *self = g_self;
  controllerJairReset(self);
}

bool controllerJairTest(controllerJair_t *self)
{
  (void)self;
  return true;
}

/* ============================================================== */
/*  Funcion principal                                             */
/* ============================================================== */
void controllerJair(controllerJair_t   *self,
                    control_t          *control,
                    const setpoint_t   *setpoint,
                    const sensorData_t *sensors,
                    const state_t      *state,
                    const uint32_t      tick)
{
  if (!RATE_DO_EXECUTE(ATTITUDE_RATE, tick)) {
    return;
  }

  const float dt = 1.0f / ATTITUDE_RATE;

  /* ============================================================
   * 0. SEGURIDAD
   * ============================================================ */
  if (setpoint->mode.z == modeDisable && setpoint->thrust < 1000) {
    control->controlMode = controlModeForceTorque;
    control->thrustSi    = 0.0f;
    control->torque[0]   = 0.0f;
    control->torque[1]   = 0.0f;
    control->torque[2]   = 0.0f;
    controllerJairReset(self);
    return;
  }

  /* ============================================================
   * 1. MEDICIONES
   * ============================================================ */
  /* Posicion y velocidad lineal */
  struct vec xi1 = mkvec(state->position.x,
                         state->position.y,
                         state->position.z);
  struct vec xi2 = mkvec(state->velocity.x,
                         state->velocity.y,
                         state->velocity.z);

  /* Orientacion: angulos de Euler desde cuaternion */
  struct quat q_cur = mkquat(state->attitudeQuaternion.x,
                             state->attitudeQuaternion.y,
                             state->attitudeQuaternion.z,
                             state->attitudeQuaternion.w);
  struct vec eta1 = quat2rpy(q_cur);
  LOG_VEC3(self, rpy, eta1);

  /* Velocidad angular del giroscopio */
  struct vec eta2 = mkvec(radians(sensors->gyro.x),
                          radians(sensors->gyro.y),
                          radians(sensors->gyro.z));

  /* ============================================================
   * 2. OBSERVADOR SM-FTO (ec.5)
   *
   * Integrador Euler hacia adelante:
   *   xi_hat1[k+1]  = xi_hat1[k]  + dt*(xi_hat2  + K1*phi1(e_xi))
   *   xi_hat2[k+1]  = xi_hat2[k]  + dt*(f_xi + xi_hat3 + K2*phi2(e_xi))
   *   xi_hat3[k+1]  = xi_hat3[k]  + dt*(K3*phi3(e_xi))
   *   eta_hat1[k+1] = eta_hat1[k] + dt*(eta_hat2 + K4*phi1(e_eta))
   *   eta_hat2[k+1] = eta_hat2[k] + dt*(f_eta + eta_hat3 + K5*phi2(e_eta))
   *   eta_hat3[k+1] = eta_hat3[k] + dt*(K6*phi3(e_eta))
   *
   * donde:
   *   f_xi  = gξ(η1)*um - G - Λξ*xi_hat2  (necesita um del ciclo actual)
   *   f_eta = J*tau + Ξ*wη(eta_hat2) - Λη*eta_hat2  (necesita tau)
   *
   * PROBLEMA: f_xi y f_eta dependen de tau y um que se calculan
   * DESPUES del observador. Para evitar la dependencia ciclica,
   * se usa la aproximacion de que la entrada del modelo (f) se
   * calcula con los valores del ciclo ANTERIOR de xi_hat2 y eta_hat2,
   * lo cual es valido para dt pequeno (1 ms).
   *
   * El observador se inicializa con las mediciones reales y
   * converge rapidamente (tiempo finito T0).
   *
   * NOTA: obs_enabled permite activar el observador de forma
   * gradual desde cfclient para evitar transitorios bruscos.
   * ============================================================ */

  /* Errores de salida del observador */
  struct vec e_xi      = vsub(xi1,  self->xi_hat1);
  struct vec e_eta_obs = vsub(eta1, self->eta_hat1);

  if (self->obs_enabled) {
    /*
     * INICIALIZACION EN FLANCO 0->1:
     * Cuando obs_enabled pasa de 0 a 1, xi_hat y eta_hat pueden
     * estar en cero (del reset) mientras las mediciones reales son
     * no-nulas. Eso produce errores de observador enormes que
     * disparan phi3=sign(e) con máxima ganancia durante muchos ciclos,
     * causando las conmutaciones observadas.
     *
     * Solución: al detectar el flanco, copiar las mediciones actuales
     * al estado del observador y resetear v_smc para evitar el
     * transitorio. El observador arranca con error cero y converge
     * suavemente desde ese punto.
     */
    if (!self->obs_prev_enabled) {
      self->xi_hat1  = xi1;
      self->xi_hat2  = xi2;
      self->xi_hat3  = vzero();
      self->eta_hat1 = eta1;
      self->eta_hat2 = eta2;
      self->eta_hat3 = vzero();
      self->v_smc    = vzero();
      /* Recalcular errores con estado recien inicializado */
      e_xi      = vzero();
      e_eta_obs = vzero();
    }
    self->obs_prev_enabled = 1;

    /* Inyecciones no lineales del observador (ec.5):
     * phi1(e) = ceil(e)^(2/3) — continua, OK sin suavizado
     * phi2(e) = ceil(e)^(1/3) — continua, OK sin suavizado
     * phi3(e) = sign(e)       — DISCONTINUA: genera chattering.
     *   Se reemplaza por sat(e/delta_obs): continuo, misma convergencia
     *   a tiempo finito para errores |e| > delta_obs.               */
    struct vec phi1_xi  = vceil_pow(e_xi,      2.0f/3.0f);
    struct vec phi2_xi  = vceil_pow(e_xi,      1.0f/3.0f);
    struct vec phi3_xi  = vceil_pow_smooth(e_xi,       0.0f, self->obs_delta_xi);
    struct vec phi1_eta = vceil_pow(e_eta_obs, 2.0f/3.0f);
    struct vec phi2_eta = vceil_pow(e_eta_obs, 1.0f/3.0f);
    struct vec phi3_eta = vceil_pow_smooth(e_eta_obs,  0.0f, self->obs_delta_eta);

    /*
     * f_xi: dinámica conocida del subsistema lineal.
     * Con Lambda_xi=0 (defecto): f_xi = 0.
     * xi_hat3 compensa gξ*um - G (perturbación efectiva).
     */
    struct vec f_xi  = vneg(veltmul(self->Lambda_xi,  self->xi_hat2));
    struct vec f_eta = vsub(coriolis_term(self->eta_hat2, self->Xi_coef),
                            veltmul(self->Lambda_eta, self->eta_hat2));

    /* Integración Euler — subsistema lineal (ec.5a-5c) */
    self->xi_hat1 = vadd(self->xi_hat1,
      vscl(dt, vadd3(self->xi_hat2,
                     f_xi,
                     veltmul(self->Kobs_1, phi1_xi))));
    self->xi_hat2 = vadd(self->xi_hat2,
      vscl(dt, vadd3(vadd(f_xi, self->xi_hat3),
                     veltmul(self->Kobs_2, phi2_xi),
                     vzero())));
    self->xi_hat3 = vadd(self->xi_hat3,
      vscl(dt, veltmul(self->Kobs_3, phi3_xi)));

    /* Integración Euler — subsistema angular (ec.5d-5f) */
    self->eta_hat1 = vadd(self->eta_hat1,
      vscl(dt, vadd(self->eta_hat2,
                    veltmul(self->Kobs_4, phi1_eta))));
    self->eta_hat2 = vadd(self->eta_hat2,
      vscl(dt, vadd3(vadd(f_eta, self->eta_hat3),
                     veltmul(self->Kobs_5, phi2_eta),
                     vzero())));
    self->eta_hat3 = vadd(self->eta_hat3,
      vscl(dt, veltmul(self->Kobs_6, phi3_eta)));

    /* Saturar estimaciones de perturbación con límites físicos:
     * xi_hat3  en [m/s²]: perturbación de aceleración lineal, max ~5 m/s²
     * eta_hat3 en [N·m]:  perturbación de torque, saturar por eje:
     *   phi/theta: 50% de tau_roll_max = 0.50 * 0.00986 = 0.00493 N·m
     *   psi:       50% de tau_yaw_max  = 0.50 * 0.00150 = 0.00075 N·m  */
    self->xi_hat3  = vclampscl(self->xi_hat3, -5.0f, 5.0f);
    self->eta_hat3.x = clampf(self->eta_hat3.x, -0.00493f, 0.00493f);
    self->eta_hat3.y = clampf(self->eta_hat3.y, -0.00493f, 0.00493f);
    self->eta_hat3.z = clampf(self->eta_hat3.z, -0.00075f, 0.00075f);

  } else {
    self->obs_prev_enabled = 0;
    /* Observador desactivado: seguir mediciones directamente     */
    self->xi_hat1  = xi1;
    self->xi_hat2  = xi2;
    self->xi_hat3  = vzero();
    self->eta_hat1 = eta1;
    self->eta_hat2 = eta2;
    self->eta_hat3 = vzero();
  }

  /* Volcar estimaciones al log */
  LOG_VEC3(self, dxi,  self->xi_hat3);
  LOG_VEC3(self, deta, self->eta_hat3);
  
    /* Velocidad y posición: Kalman vs observador (para guardar_experimento.py) */
  self->log_vel_kx = xi2.x;            /* Kalman directo */
  self->log_vel_ky = xi2.y;
  self->log_vel_kz = xi2.z;
  self->log_vel_ox = self->xi_hat2.x;  /* observador SM-FTO */
  self->log_vel_oy = self->xi_hat2.y;
  self->log_vel_oz = self->xi_hat2.z;
  self->log_pos_kx = xi1.x;            /* Kalman directo */
  self->log_pos_ky = xi1.y;
  self->log_pos_kz = xi1.z;
  self->log_pos_ox = self->xi_hat1.x;  /* observador SM-FTO */
  self->log_pos_oy = self->xi_hat1.y;
  self->log_pos_oz = self->xi_hat1.z;


  /* ============================================================
   * 3. GUINADA DESEADA psi_d
   * ============================================================ */
  float psi_d = 0.0f;
  if (setpoint->mode.yaw == modeVelocity) {
    psi_d = eta1.z + radians(setpoint->attitudeRate.yaw) * dt;
    psi_d = wrap_angle(psi_d);
  } else if (setpoint->mode.yaw == modeAbs) {
    psi_d = radians(setpoint->attitude.yaw);
  } else if (setpoint->mode.quat == modeAbs) {
    struct quat q_sp = mkquat(setpoint->attitudeQuaternion.x,
                              setpoint->attitudeQuaternion.y,
                              setpoint->attitudeQuaternion.z,
                              setpoint->attitudeQuaternion.w);
    psi_d = quat2rpy(q_sp).z;
  }

  /* ============================================================
   * 4. CONTROL DE POSICION (ec.12 + 11a)
   *
   * nu_bar = k0xi*int(exi) + k1xi*exi + k2xi*eps_xi
   * nu = nu_bar + Lambda_xi*xi2 + xi_dd - xi_hat3  (ec.11a)
   *
   * Las estimaciones del observador se usan en xi_hat2 (velocidad
   * estimada, mejor que la medida directa con ruido) y xi_hat3
   * (compensacion de perturbacion).
   * ============================================================ */
  float phi_star   = 0.0f;
  float theta_star = 0.0f;

  if (   setpoint->mode.x == modeAbs
      || setpoint->mode.y == modeAbs
      || setpoint->mode.z == modeAbs) {

    /*
     * GUARDIAN DE ESTIMADOR:
     * Si el estimador de posicion no ha convergido (Z fuera de rango
     * fisico razonable, p.ej. > 100 m), tratar como modo manual para
     * evitar que el controlador de posicion genere comandos absurdos.
     * Esto ocurre cuando se usa "Position hc" en cfclient sin deck de
     * posicionamiento: el filtro Kalman reporta valores basura.
     * En ese caso se ignora el setpoint de posicion y se usa la
     * referencia de actitud del setpoint de forma plana.
     */
    if (fabsf(xi1.z) > 100.0f || fabsf(xi1.x) > 100.0f || fabsf(xi1.y) > 100.0f) {
      float thrust_norm  = (float)setpoint->thrust / 65535.0f;
      control->thrustSi  = thrust_norm * CF2_MAX_THRUST_N;
      phi_star   =  radians(setpoint->attitude.roll);
      theta_star = -radians(setpoint->attitude.pitch);
      self->log_thrustSi   = control->thrustSi;
      self->log_phi_star   = phi_star;
      self->log_theta_star = theta_star;
      /* Saltar directamente al control de orientacion */
      goto attitude_control;
    }

    struct vec xi_d  = mkvec(setpoint->position.x,
                             setpoint->position.y,
                             setpoint->position.z);

    /* Filtro anti-spike del estimador Kalman:
     * Si la posicion reportada salta mas de 50cm en un ciclo
     * (imposible fisicamente a 1kHz), ignorar esa medicion.
     * Esto previo el crash del segundo vuelo donde ex=1.14m
     * aparecio en un solo ciclo por ruido del estimador.     */
    if (self->prev_valid) {
      struct vec jump = vsub(xi1, self->xi_hat1);
      float jump_mag = sqrtf(jump.x*jump.x + jump.y*jump.y + jump.z*jump.z);
      if (jump_mag > 0.50f) {
        /* Spike detectado: usar ultima posicion conocida    */
        controllerJairReset(self);
        control->controlMode = controlModeForceTorque;
        control->thrustSi    = self->mass * GRAVITY_MAGNITUDE;
        control->torque[0]   = 0.0f;
        control->torque[1]   = 0.0f;
        control->torque[2]   = 0.0f;
        return;
      }
    }
    struct vec xi_d2 = mkvec(setpoint->velocity.x,
                             setpoint->velocity.y,
                             setpoint->velocity.z);

    /* SIEMPRE usar velocidad del estimador Kalman para eps_xi.
     * xi_hat2 del observador SM-FTO tiene el mismo problema que
     * eta_hat2: phi2 amplifica ruido. xi_hat3 sí se usa en nu. */
    struct vec xi2_est = xi2;  /* siempre Kalman directo */

    struct vec e_xi_ctrl = vclampscl(vsub(xi_d,  xi1),
                                     -self->Kpos_P_limit,
                                      self->Kpos_P_limit);
    struct vec eps_xi    = vclampscl(vsub(xi_d2, xi2_est),
                                     -self->Kpos_D_limit,
                                      self->Kpos_D_limit);
    LOG_VEC3(self, epos, e_xi_ctrl);
    LOG_VEC3(self, evel, eps_xi);

    self->i_error_pos = vadd(self->i_error_pos, vscl(dt, e_xi_ctrl));
    self->i_error_pos = vclampscl(self->i_error_pos,
                                  -self->Kpos_I_limit,
                                   self->Kpos_I_limit);

    /* nu_bar = k1xi*exi + k2xi*eps_xi + k0xi*int(exi)  (ec.12) */
    struct vec nu_bar = vadd3(veltmul(self->Kpos_P, e_xi_ctrl),
                              veltmul(self->Kpos_D, eps_xi),
                              veltmul(self->Kpos_I, self->i_error_pos));

    /* Aceleracion de referencia xi_dd */
    struct vec xi_dd;
    {
      bool sp_has_acc =
        (fabsf(setpoint->acceleration.x) > 1e-4f) ||
        (fabsf(setpoint->acceleration.y) > 1e-4f) ||
        (fabsf(setpoint->acceleration.z) > 1e-4f);
      if (sp_has_acc) {
        xi_dd = mkvec(setpoint->acceleration.x,
                      setpoint->acceleration.y,
                      setpoint->acceleration.z);
      } else if (self->prev_valid) {
        xi_dd = vclampscl(
          vscl(1.0f / dt, vsub(xi_d2, self->prev_xi_d2)),
          -5.0f, 5.0f);
      } else {
        xi_dd = vzero();
      }
    }
    self->prev_xi_d2 = xi_d2;
    LOG_VEC3(self, xidd, xi_dd);

    /* nu = nu_bar + Lambda_xi*xi2 + xi_dd - xi_hat3  (ec.11a)  */
    struct vec nu = vsub(
      vadd3(nu_bar, veltmul(self->Lambda_xi, xi2_est), xi_dd),
      self->xi_hat3);
    LOG_VEC3(self, nu, nu);

    /* Referencias angulares (ec.10) */
    float nu_z_g = nu.z + GRAVITY_MAGNITUDE;
    float u_m    = sqrtf(nu.x*nu.x + nu.y*nu.y + nu_z_g*nu_z_g);

    if (u_m > 1e-6f) {
      float sp  = sinf(psi_d);
      float cp  = cosf(psi_d);
      float arg = clampf((nu.x*sp - nu.y*cp) / u_m, -1.0f, 1.0f);
      phi_star   = asinf(arg);
      theta_star = atanf((nu.x*cp + nu.y*sp) / nu_z_g);
      /* Limite de angulo anti-volteo: ajustable desde cfclient
       * como ctrlJair.max_tilt. Defecto 20 deg (0.3491 rad).    */
      phi_star   = clampf(phi_star,   -self->max_tilt_rad, self->max_tilt_rad);
      theta_star = clampf(theta_star, -self->max_tilt_rad, self->max_tilt_rad);
    }

    /* Thrust con feedforward gravitacional:
     * thrustSi = mass * u_m donde u_m incluye g via nu_z+g.
     * Si la masa configurada < masa real, el quad no sube.
     * Verificar: hovering debe ocurrir con thrustSi ~ mass*g. */
    control->thrustSi = self->mass * u_m;

    /* Saturacion de seguridad: nunca mas de 4x el peso del dron.
     * Evita que un estimador desbocado sature los motores.        */
    /* CF2_MAX_THRUST_N = 0.480 N (ya definido arriba) */
    if (control->thrustSi > CF2_MAX_THRUST_N) {
      control->thrustSi = CF2_MAX_THRUST_N;
    }

    if (control->thrustSi < 0.01f) {
      controllerJairReset(self);
    }

  } else {
    /* Modo manual */
    float thrust_norm  = (float)setpoint->thrust / 65535.0f;
    control->thrustSi  = thrust_norm * CF2_MAX_THRUST_N;
    phi_star   =  radians(setpoint->attitude.roll);
    theta_star = -radians(setpoint->attitude.pitch);
    self->log_nu_x   = 0.0f; self->log_nu_y   = 0.0f; self->log_nu_z   = 0.0f;
    self->log_xidd_x = 0.0f; self->log_xidd_y = 0.0f; self->log_xidd_z = 0.0f;
  }

  self->log_thrustSi   = control->thrustSi;
  self->log_phi_star   = phi_star;
  self->log_theta_star = theta_star;

  /* Etiqueta de salto para el guardian de estimador */
  attitude_control:;

  /* ============================================================
   * 5. CONTROL DE ORIENTACION STSMC (ec.18 + 11b)
   *
   * s = eps_eta + k0 * ceil(e_eta)^(2/3)
   * tau_bar = v - k1 * ceil(s)^(1/2)
   * v_dot   = -k2 * sign(s)   [suavizado con delta_s]
   *
   * tau (ec.11b) = tau_bar - Xi*wη(eta2) + Lambda_eta*eta2
   *               + eta_dd - eta_hat3
   *
   * NOTA: el paper escribe tau = J^-1*(tau_bar + ...).
   * Aqui se omite J^-1 porque control->torque[] espera [N*m]
   * (ver nota al inicio del archivo).
   * ============================================================ */

  /* Referencias angulares eta1_d = (phi*, theta*, psi_d) */
  struct vec eta1_d = mkvec(phi_star, theta_star, psi_d);
  LOG_VEC3(self, rpyd, eta1_d);

  /*
   * eta_dot_d por filtro IIR de primer orden.
   * Se limita la derivada cruda a +-10 rad/s antes del filtro
   * para evitar spikes cuando phi_star/theta_star cambian bruscamente.
   */
  if (self->prev_valid) {
    struct vec eta1d_diff = vsub(eta1_d, self->prev_eta1_d);
    /* Wrap del componente yaw para diferencias angulares correctas */
    eta1d_diff.z = wrap_angle(eta1d_diff.z);
    struct vec eta2d_raw = vclampscl(vscl(1.0f/dt, eta1d_diff),
                                     -10.0f, 10.0f);
    float a = self->eta2d_filter_alpha;
    self->eta2d_filt = vadd(vscl(a, self->eta2d_filt),
                            vscl(1.0f - a, eta2d_raw));
  } else {
    /* PRIMER CICLO VALIDO: inicializar eta2d_filt con el giroscopio
     * actual para que eps_eta = eta2 - eta2d_filt = 0 en el arranque.
     * Sin esto, la vibracion mecanica del encendido de motores
     * (tipicamente 0.5-2 rad/s en el gyro) genera eps_eta != 0 y
     * el STSMC aplica torque asimetrico desde el ciclo 1, desplazando
     * el quad. El PID de Bitcraze no tiene este problema porque usa
     * error de angulo (rpy_real - rpy_sp ~ 0 en tierra), no gyro directo.
     * Inicializar aqui garantiza eps_eta=0 y s=k0*e^(2/3)~0 al arranque. */
    self->eta2d_filt = eta2;
  }
  self->prev_eta1_d = eta1_d;
  self->prev_valid  = true;

  /* SIEMPRE usar el giroscopio directo para eps_eta.
   * eta_hat2 NO se usa en el controlador: phi2(e)=e^(1/3) amplifica
   * el ruido IMU (~0.5°) a ±0.20 rad/s en eta_hat2, lo que hace
   * que eps_eta fluctúe fuera de delta_s y cause chattering.
   * eta_hat3 (perturbación estimada) SÍ se usa en tau_num. */
  struct vec eta2_ctrl = eta2;  /* siempre giroscopio directo */

  /*
   * Error de orientacion con wrap en yaw (evita giro continuo).
   */
  struct vec e_eta;
  e_eta.x = wrap_angle(eta1.x - eta1_d.x);
  e_eta.y = wrap_angle(eta1.y - eta1_d.y);
  e_eta.z = wrap_angle(eta1.z - eta1_d.z);

  struct vec eps_eta = vsub(eta2_ctrl, self->eta2d_filt);
  LOG_VEC3(self, eatt,  e_eta);
  LOG_VEC3(self, erate, eps_eta);

  /*
   * Superficie: s = eps_eta + k0 * ceil(e_eta)^(2/3)
   * Se suaviza ceil(e_eta)^(2/3) con delta_e.
   */
  struct vec s_smc = vadd(eps_eta,
                          veltmul(self->Ksmc_k0,
                                  vceil_pow_smooth(e_eta, 2.0f/3.0f,
                                                   self->Ksmc_delta_e)));
  LOG_VEC3(self, s, s_smc);

  /*
   * Integrador v: v[k] = v[k-1] - k2 * sat(s/delta_s) * dt
   */
  self->v_smc = vsub(self->v_smc,
                     vscl(dt, veltmul(self->Ksmc_k2,
                                      vceil_pow_smooth(s_smc, 0.0f,
                                                       self->Ksmc_delta_s))));
  /* Saturacion anti-windup por eje: yaw tiene limite mucho menor  */
  self->v_smc = mkvec(
    clampf(self->v_smc.x, -self->Ksmc_max_v.x, self->Ksmc_max_v.x),
    clampf(self->v_smc.y, -self->Ksmc_max_v.y, self->Ksmc_max_v.y),
    clampf(self->v_smc.z, -self->Ksmc_max_v.z, self->Ksmc_max_v.z));
  LOG_VEC3(self, v, self->v_smc);

  /*
   * tau_bar = v - k1 * ceil(s)^(1/2)
   */
  struct vec tau_bar = vsub(self->v_smc,
                            veltmul(self->Ksmc_k1,
                                    vceil_pow_smooth(s_smc, 0.5f,
                                                     self->Ksmc_delta_s)));

  /*
   * tau (ec.11b):
   *   tau = tau_bar - Xi*wη(eta2) + Lambda_eta*eta2
   *         + eta_dd - eta_hat3
   *
   * - Xi*wη: compensacion de Coriolis (desactivado por defecto)
   * - Lambda_eta*eta2: amortiguacion (desactivada por defecto)
   * - eta_hat3: compensacion de perturbacion angular del observador
   * - eta_dd: se omite (aproximado como cero, el STSMC lo compensa)
   */
  struct vec xi_w    = coriolis_term(eta2, self->Xi_coef);
  struct vec tau_num = vsub(
    vadd3(tau_bar,
          vneg(xi_w),
          veltmul(self->Lambda_eta, eta2)),
    self->eta_hat3);

  /* Saturación por eje con límites físicos reales */
  struct vec tau = mkvec(
    clampf(tau_num.x, -MAX_TORQUE_ROLL_NM, MAX_TORQUE_ROLL_NM),
    clampf(tau_num.y, -MAX_TORQUE_ROLL_NM, MAX_TORQUE_ROLL_NM),
    clampf(tau_num.z, -MAX_TORQUE_YAW_NM,  MAX_TORQUE_YAW_NM));
  LOG_VEC3(self, tau, tau);

  /* ============================================================
   * 6. SALIDAS AL FIRMWARE
   * ============================================================ */
  control->controlMode = controlModeForceTorque;
  control->torque[0]   = tau.x;
  control->torque[1]   = tau.y;
  control->torque[2]   = tau.z;
}

/* ============================================================== */
/*  WRAPPER FIRMWARE CRAZYFLIE                                    */
/* ============================================================== */
#ifdef CRAZYFLIE_FW

#include "param.h"
#include "log.h"

void controllerJairFirmwareInit(void)
{
  controllerJairInit(&g_self);
}

bool controllerJairFirmwareTest(void)
{
  return controllerJairTest(&g_self);
}

void controllerJairFirmware(control_t          *control,
                            const setpoint_t   *setpoint,
                            const sensorData_t *sensors,
                            const state_t      *state,
                            const uint32_t      tick)
{
  controllerJair(&g_self, control, setpoint, sensors, state, tick);
}

/* -------------------------------------------------------------- */
/*  PARAMETROS (ajustables desde cfclient / cflib en vuelo)       */
/* -------------------------------------------------------------- */
PARAM_GROUP_START(ctrlJair)
/* Fisicos */
PARAM_ADD(PARAM_FLOAT, mass,        &g_self.mass)
PARAM_ADD(PARAM_FLOAT, Jx,          &g_self.J.x)
PARAM_ADD(PARAM_FLOAT, Jy,          &g_self.J.y)
PARAM_ADD(PARAM_FLOAT, Jz,          &g_self.J.z)
/* Posicion PID */
PARAM_ADD(PARAM_FLOAT, Kpos_Px,     &g_self.Kpos_P.x)
PARAM_ADD(PARAM_FLOAT, Kpos_Py,     &g_self.Kpos_P.y)
PARAM_ADD(PARAM_FLOAT, Kpos_Pz,     &g_self.Kpos_P.z)
PARAM_ADD(PARAM_FLOAT, max_tilt,    &g_self.max_tilt_rad)
PARAM_ADD(PARAM_FLOAT, Kpos_Plim,   &g_self.Kpos_P_limit)
PARAM_ADD(PARAM_FLOAT, Kpos_Dx,     &g_self.Kpos_D.x)
PARAM_ADD(PARAM_FLOAT, Kpos_Dy,     &g_self.Kpos_D.y)
PARAM_ADD(PARAM_FLOAT, Kpos_Dz,     &g_self.Kpos_D.z)
PARAM_ADD(PARAM_FLOAT, Kpos_Dlim,   &g_self.Kpos_D_limit)
PARAM_ADD(PARAM_FLOAT, Kpos_Ix,     &g_self.Kpos_I.x)
PARAM_ADD(PARAM_FLOAT, Kpos_Iy,     &g_self.Kpos_I.y)
PARAM_ADD(PARAM_FLOAT, Kpos_Iz,     &g_self.Kpos_I.z)
PARAM_ADD(PARAM_FLOAT, Kpos_Ilim,   &g_self.Kpos_I_limit)
/* STSMC orientacion */
PARAM_ADD(PARAM_FLOAT, k0_phi,      &g_self.Ksmc_k0.x)
PARAM_ADD(PARAM_FLOAT, k0_tht,      &g_self.Ksmc_k0.y)
PARAM_ADD(PARAM_FLOAT, k0_psi,      &g_self.Ksmc_k0.z)
PARAM_ADD(PARAM_FLOAT, k1_phi,      &g_self.Ksmc_k1.x)
PARAM_ADD(PARAM_FLOAT, k1_tht,      &g_self.Ksmc_k1.y)
PARAM_ADD(PARAM_FLOAT, k1_psi,      &g_self.Ksmc_k1.z)
PARAM_ADD(PARAM_FLOAT, k2_phi,      &g_self.Ksmc_k2.x)
PARAM_ADD(PARAM_FLOAT, k2_tht,      &g_self.Ksmc_k2.y)
PARAM_ADD(PARAM_FLOAT, k2_psi,      &g_self.Ksmc_k2.z)
PARAM_ADD(PARAM_FLOAT, maxv_phi,    &g_self.Ksmc_max_v.x)
PARAM_ADD(PARAM_FLOAT, maxv_tht,    &g_self.Ksmc_max_v.y)
PARAM_ADD(PARAM_FLOAT, maxv_psi,    &g_self.Ksmc_max_v.z)
PARAM_ADD(PARAM_FLOAT, delta_s,     &g_self.Ksmc_delta_s)
PARAM_ADD(PARAM_FLOAT, delta_e,     &g_self.Ksmc_delta_e)
/* Observador */
PARAM_ADD(PARAM_UINT8, obs_en,      &g_self.obs_enabled)
PARAM_ADD(PARAM_FLOAT, obs_dxi,     &g_self.obs_delta_xi)
PARAM_ADD(PARAM_FLOAT, obs_deta,    &g_self.obs_delta_eta)
PARAM_ADD(PARAM_FLOAT, K1x,         &g_self.Kobs_1.x)
PARAM_ADD(PARAM_FLOAT, K1y,         &g_self.Kobs_1.y)
PARAM_ADD(PARAM_FLOAT, K1z,         &g_self.Kobs_1.z)
PARAM_ADD(PARAM_FLOAT, K2x,         &g_self.Kobs_2.x)
PARAM_ADD(PARAM_FLOAT, K2y,         &g_self.Kobs_2.y)
PARAM_ADD(PARAM_FLOAT, K2z,         &g_self.Kobs_2.z)
PARAM_ADD(PARAM_FLOAT, K3x,         &g_self.Kobs_3.x)
PARAM_ADD(PARAM_FLOAT, K3y,         &g_self.Kobs_3.y)
PARAM_ADD(PARAM_FLOAT, K3z,         &g_self.Kobs_3.z)
PARAM_ADD(PARAM_FLOAT, K4x,         &g_self.Kobs_4.x)
PARAM_ADD(PARAM_FLOAT, K4y,         &g_self.Kobs_4.y)
PARAM_ADD(PARAM_FLOAT, K4z,         &g_self.Kobs_4.z)
PARAM_ADD(PARAM_FLOAT, K5x,         &g_self.Kobs_5.x)
PARAM_ADD(PARAM_FLOAT, K5y,         &g_self.Kobs_5.y)
PARAM_ADD(PARAM_FLOAT, K5z,         &g_self.Kobs_5.z)
PARAM_ADD(PARAM_FLOAT, K6x,         &g_self.Kobs_6.x)
PARAM_ADD(PARAM_FLOAT, K6y,         &g_self.Kobs_6.y)
PARAM_ADD(PARAM_FLOAT, K6z,         &g_self.Kobs_6.z)
/* Filtro eta_dot_d */
PARAM_ADD(PARAM_FLOAT, eta2d_alfa,  &g_self.eta2d_filter_alpha)
/* Amortiguacion (desactivada por defecto) */
PARAM_ADD(PARAM_FLOAT, Lxi_x,       &g_self.Lambda_xi.x)
PARAM_ADD(PARAM_FLOAT, Lxi_y,       &g_self.Lambda_xi.y)
PARAM_ADD(PARAM_FLOAT, Lxi_z,       &g_self.Lambda_xi.z)
PARAM_ADD(PARAM_FLOAT, Leta_x,      &g_self.Lambda_eta.x)
PARAM_ADD(PARAM_FLOAT, Leta_y,      &g_self.Lambda_eta.y)
PARAM_ADD(PARAM_FLOAT, Leta_z,      &g_self.Lambda_eta.z)
PARAM_GROUP_STOP(ctrlJair)

/* -------------------------------------------------------------- */
/*  LOG                                                           */
/*  Todos los ADDRESS son &g_self.log_campo (float escalar)       */
/* -------------------------------------------------------------- */
LOG_GROUP_START(ctrlJair)
LOG_ADD(LOG_FLOAT, thrustSi,    &g_self.log_thrustSi)
LOG_ADD(LOG_FLOAT, phi_star,    &g_self.log_phi_star)
LOG_ADD(LOG_FLOAT, theta_star,  &g_self.log_theta_star)
LOG_ADD(LOG_FLOAT, rpy_x,       &g_self.log_rpy_x)
LOG_ADD(LOG_FLOAT, rpy_y,       &g_self.log_rpy_y)
LOG_ADD(LOG_FLOAT, rpy_z,       &g_self.log_rpy_z)
LOG_ADD(LOG_FLOAT, rpyd_x,      &g_self.log_rpyd_x)
LOG_ADD(LOG_FLOAT, rpyd_y,      &g_self.log_rpyd_y)
LOG_ADD(LOG_FLOAT, rpyd_z,      &g_self.log_rpyd_z)
LOG_ADD(LOG_FLOAT, epos_x,      &g_self.log_epos_x)
LOG_ADD(LOG_FLOAT, epos_y,      &g_self.log_epos_y)
LOG_ADD(LOG_FLOAT, epos_z,      &g_self.log_epos_z)
LOG_ADD(LOG_FLOAT, evel_x,      &g_self.log_evel_x)
LOG_ADD(LOG_FLOAT, evel_y,      &g_self.log_evel_y)
LOG_ADD(LOG_FLOAT, evel_z,      &g_self.log_evel_z)
LOG_ADD(LOG_FLOAT, eatt_x,      &g_self.log_eatt_x)
LOG_ADD(LOG_FLOAT, eatt_y,      &g_self.log_eatt_y)
LOG_ADD(LOG_FLOAT, eatt_z,      &g_self.log_eatt_z)
LOG_ADD(LOG_FLOAT, erate_x,     &g_self.log_erate_x)
LOG_ADD(LOG_FLOAT, erate_y,     &g_self.log_erate_y)
LOG_ADD(LOG_FLOAT, erate_z,     &g_self.log_erate_z)
LOG_ADD(LOG_FLOAT, s_x,         &g_self.log_s_x)
LOG_ADD(LOG_FLOAT, s_y,         &g_self.log_s_y)
LOG_ADD(LOG_FLOAT, s_z,         &g_self.log_s_z)
LOG_ADD(LOG_FLOAT, v_x,         &g_self.log_v_x)
LOG_ADD(LOG_FLOAT, v_y,         &g_self.log_v_y)
LOG_ADD(LOG_FLOAT, v_z,         &g_self.log_v_z)
LOG_ADD(LOG_FLOAT, tau_x,       &g_self.log_tau_x)
LOG_ADD(LOG_FLOAT, tau_y,       &g_self.log_tau_y)
LOG_ADD(LOG_FLOAT, tau_z,       &g_self.log_tau_z)
LOG_ADD(LOG_FLOAT, nu_x,        &g_self.log_nu_x)
LOG_ADD(LOG_FLOAT, nu_y,        &g_self.log_nu_y)
LOG_ADD(LOG_FLOAT, nu_z,        &g_self.log_nu_z)
LOG_ADD(LOG_FLOAT, xidd_x,      &g_self.log_xidd_x)
LOG_ADD(LOG_FLOAT, xidd_y,      &g_self.log_xidd_y)
LOG_ADD(LOG_FLOAT, xidd_z,      &g_self.log_xidd_z)
LOG_ADD(LOG_FLOAT, dxi_x,       &g_self.log_dxi_x)
LOG_ADD(LOG_FLOAT, dxi_y,       &g_self.log_dxi_y)
LOG_ADD(LOG_FLOAT, dxi_z,       &g_self.log_dxi_z)
LOG_ADD(LOG_FLOAT, deta_x,      &g_self.log_deta_x)
LOG_ADD(LOG_FLOAT, deta_y,      &g_self.log_deta_y)
LOG_ADD(LOG_FLOAT, deta_z,      &g_self.log_deta_z)
LOG_ADD(LOG_FLOAT, vel_kx,      &g_self.log_vel_kx)
LOG_ADD(LOG_FLOAT, vel_ky,      &g_self.log_vel_ky)
LOG_ADD(LOG_FLOAT, vel_kz,      &g_self.log_vel_kz)
LOG_ADD(LOG_FLOAT, vel_ox,      &g_self.log_vel_ox)
LOG_ADD(LOG_FLOAT, vel_oy,      &g_self.log_vel_oy)
LOG_ADD(LOG_FLOAT, vel_oz,      &g_self.log_vel_oz)
LOG_ADD(LOG_FLOAT, pos_kx,      &g_self.log_pos_kx)
LOG_ADD(LOG_FLOAT, pos_ky,      &g_self.log_pos_ky)
LOG_ADD(LOG_FLOAT, pos_kz,      &g_self.log_pos_kz)
LOG_ADD(LOG_FLOAT, pos_ox,      &g_self.log_pos_ox)
LOG_ADD(LOG_FLOAT, pos_oy,      &g_self.log_pos_oy)
LOG_ADD(LOG_FLOAT, pos_oz,      &g_self.log_pos_oz)
LOG_GROUP_STOP(ctrlJair)

#endif /* CRAZYFLIE_FW */
