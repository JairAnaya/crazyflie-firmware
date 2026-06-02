/*
 * controller_jair.h
 *
 * Implementacion fiel al articulo:
 *   J. Anaya, A.E. Dzul, H. Rios, S. Encina-Espino
 *   "Estrategia de control basada en Modos Deslizantes para el
 *    vuelo en formacion de Quad-Rotors", COMRob 2025
 *
 * Arquitectura:
 *   POSICION  : PID (ec.12) + control virtual nu (ec.11a)
 *   ORIENTACION: STSMC (ec.18) con tau (ec.11b)
 *   OBSERVADOR : SM-FTO (ec.5) para xi y eta
 *
 * Nota sobre unidades:
 *   Las ganancias del observador y del STSMC en el paper estan en
 *   escala de aceleracion [rad/s²] porque la ec.11b incluye J^-1.
 *   En el firmware Crazyflie, control->torque[] espera [N·m] y
 *   power_distribution aplica J^-1 internamente. Por tanto:
 *     - tau_bar se entrega en [N·m] directamente
 *     - Las ganancias k1, k2 del STSMC se escalan por J:
 *         k1_Nm = k1_paper * J,  k2_Nm = k2_paper * J
 *     - k0 es adimensional (superficie en [rad/s])
 *
 * Campos de log declarados como float escalares (no struct vec)
 * para garantizar que LOG_ADD resuelva &g_self.campo como constante
 * en tiempo de compilacion (ARM-GCC, Crazyflie firmware).
 */

#pragma once

#include <stdbool.h>
#include <stdint.h>

#include "math3d.h"
#include "stabilizer_types.h"

typedef struct controllerJair_s {

  /* --- Parametros fisicos --------------------------------------- */
  float mass;
  struct vec J;           /* momentos de inercia [kg*m^2]          */

  /* --- Ganancias posicion PID (ec.12, Tabla IV-VI) ------------- */
  /* nu_bar = k0xi*int(exi) + k1xi*exi + k2xi*epsi_xi             */
  struct vec Kpos_P;      /* k1xi [m/s^2 / m]                      */
  struct vec Kpos_D;      /* k2xi [m/s^2 / (m/s)]                  */
  struct vec Kpos_I;      /* k0xi [m/s^2 / (m*s)]                  */
  float      Kpos_P_limit;
  float      Kpos_D_limit;
  float      Kpos_I_limit;
  float      max_tilt_rad; /* limite angulo phi* theta* [rad]   */
  struct vec Lambda_xi;   /* amortiguacion lineal (0 por defecto)   */

  /* --- Ganancias STSMC orientacion (ec.18, Tabla IV-VI) -------- */
  /* s = eps_eta + k0 * ceil(e_eta)^(2/3)                         */
  /* tau_bar = v - k1 * ceil(s)^(1/2)                             */
  /* v_dot   = -k2 * sign(s)                                      */
  /* UNIDADES en [N*m]: k1_Nm = k1_paper * J, k2_Nm = k2_paper*J */
  struct vec Ksmc_k0;     /* adimensional [rad/s / rad^(2/3)]      */
  struct vec Ksmc_k1;     /* [N*m / (rad/s)^0.5]                   */
  struct vec Ksmc_k2;     /* [N*m / s]                             */
  struct vec Ksmc_max_v;  /* saturacion anti-windup de v [N*m], por eje */
  /* Capa limite anti-chattering (sat en lugar de sign)            */
  float      Ksmc_delta_s;/* sobre s [rad/s]                       */
  float      Ksmc_delta_e;/* sobre e_eta [rad]                     */
  struct vec Lambda_eta;  /* amortiguacion angular (0 por defecto)  */
  struct vec Xi_coef;     /* Coriolis escalado [adim]               */

  /* --- Ganancias observador SM-FTO (ec.5, del paper Sec.V) ----- */
  /* phi1(e) = ceil(e)^(2/3),  phi2(e) = ceil(e)^(1/3)            */
  /* phi3(e) = sign(e)                                             */
  /* xih_dot1 = xih2 + K1*phi1(e_xi)                              */
  /* xih_dot2 = f_xi + xih3 + K2*phi2(e_xi)                       */
  /* xih_dot3 = K3*phi3(e_xi)                                      */
  /* etah_dot1 = etah2 + K4*phi1(e_eta)                           */
  /* etah_dot2 = f_eta + etah3 + K5*phi2(e_eta)                   */
  /* etah_dot3 = K6*phi3(e_eta)                                    */
  struct vec Kobs_1;      /* K1: diag(1.0766,1.0766,1.2324)        */
  struct vec Kobs_2;      /* K2: diag(0.5925,0.5925,0.7256)        */
  /* K3 y K6 se calculan como K3=D_xi*obs_margin, K6=D_eta*obs_margin
   * El usuario configura D_xi, D_eta desde cfclient según la
   * perturbación esperada. K3 y K6 se actualizan automáticamente. */
  float      obs_D_xi;    /* cota superior |d_xi_dot|  [m/s³]      */
  float      obs_D_eta;   /* cota superior |d_eta_dot| [N·m/s]     */
  float      obs_margin;  /* margen sobre D: K=D*margin (def 1.2)  */
  struct vec Kobs_3;      /* K3 = D_xi  * margin (calculado)       */
  struct vec Kobs_4;      /* K4: diag(1.4930,1.6083,1.6602)        */
  struct vec Kobs_5;      /* K5: diag(0.9675,1.0817,1.1345)        */
  struct vec Kobs_6;      /* K6 = D_eta * margin (calculado)       */
  uint8_t    obs_enabled; /* activar/desactivar observador          */
  float      obs_delta_xi;  /* capa límite phi3 observador pos [m]  */
  float      obs_delta_eta; /* capa límite phi3 observador ang [rad] */

  /* --- Estado interno del observador SM-FTO -------------------- */
  struct vec xi_hat1;     /* posicion estimada [m]                 */
  struct vec xi_hat2;     /* velocidad estimada [m/s]              */
  struct vec xi_hat3;     /* perturbacion lineal estimada [m/s^2]  */
  struct vec eta_hat1;    /* angulo estimado [rad]                 */
  struct vec eta_hat2;    /* vel. angular estimada [rad/s]         */
  struct vec eta_hat3;    /* perturbacion angular estimada [N*m]   */

  /* --- Estado interno del controlador -------------------------- */
  struct vec i_error_pos; /* integral error posicion               */
  struct vec v_smc;       /* integrador STSMC [N*m]               */
  struct vec prev_eta1_d; /* eta1_d ciclo anterior [rad]           */
  struct vec eta2d_filt;  /* eta_dot_d filtrado [rad/s]            */
  struct vec prev_xi_d2;  /* xi_d2 ciclo anterior [m/s]           */
  bool       prev_valid;
  float      eta2d_filter_alpha;
  uint8_t    obs_prev_enabled;   /* para detectar flanco 0->1      */

  /* --- Campos de log (float escalar por componente) ------------ */
  /* CRITICO: deben ser float directos para que LOG_ADD compile    */
  float log_thrustSi;
  float log_phi_star;
  float log_theta_star;
  float log_rpy_x;
  float log_rpy_y;
  float log_rpy_z;
  float log_rpyd_x;
  float log_rpyd_y;
  float log_rpyd_z;
  float log_epos_x;
  float log_epos_y;
  float log_epos_z;
  float log_evel_x;
  float log_evel_y;
  float log_evel_z;
  float log_eatt_x;
  float log_eatt_y;
  float log_eatt_z;
  float log_erate_x;
  float log_erate_y;
  float log_erate_z;
  float log_s_x;
  float log_s_y;
  float log_s_z;
  float log_v_x;
  float log_v_y;
  float log_v_z;
  float log_tau_x;
  float log_tau_y;
  float log_tau_z;
  float log_nu_x;
  float log_nu_y;
  float log_nu_z;
  float log_xidd_x;
  float log_xidd_y;
  float log_xidd_z;
  /* Perturbaciones estimadas por el observador */
  float log_dxi_x;
  float log_dxi_y;
  float log_dxi_z;
  float log_deta_x;
  float log_deta_y;
  float log_deta_z;
  /* Velocidad Kalman (real) vs observador (estimada) — para validación */
  float log_vel_kx;   /* state->velocity.x  [m/s] */
  float log_vel_ky;
  float log_vel_kz;
  float log_vel_ox;   /* xi_hat2.x observador [m/s] */
  float log_vel_oy;
  float log_vel_oz;
  /* Posición Kalman vs observador — para validación de xi_hat1 */
  float log_pos_kx;   /* state->position.x  [m] */
  float log_pos_ky;
  float log_pos_kz;
  float log_pos_ox;   /* xi_hat1.x observador [m] */
  float log_pos_oy;
  float log_pos_oz;

} controllerJair_t;

/* --- API publica --------------------------------------------- */
void controllerJairInit(controllerJair_t *self);
void controllerJairReset(controllerJair_t *self);
bool controllerJairTest(controllerJair_t *self);

void controllerJair(controllerJair_t   *self,
                    control_t          *control,
                    const setpoint_t   *setpoint,
                    const sensorData_t *sensors,
                    const state_t      *state,
                    const uint32_t      tick);

#ifdef CRAZYFLIE_FW
void controllerJairFirmwareInit(void);
bool controllerJairFirmwareTest(void);
void controllerJairFirmware(control_t          *control,
                            const setpoint_t   *setpoint,
                            const sensorData_t *sensors,
                            const state_t      *state,
                            const uint32_t      tick);
#endif
