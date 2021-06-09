#ifndef TRIGONOMETRY_H
#define TRIGONOMETRY_H
/* Include Files */

#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
/* structure type defination */
typedef struct {
	uint32_t s_buff[2];
} int64m_t;
typedef struct {
	uint32_t s_buff[3];
} int96m_t;
typedef struct {
	int32_t re;
	int32_t im;
} cint32_T;

/* State_type */
typedef enum
{
	cordic_acosine_state,
	cordic_cosine_state,
	cordic_sincos_state,
	cordic_angle_state,
	cordic_asine_state,
	cordic_atan_state,
	cordic_cexp_state,
	cordic_sine_state,
	cordic_last_state
} cordicstate;

/* Event_type */
typedef enum
{
	cordic_acosine_event,
	cordic_cosine_event,
	cordic_sincos_event,
	cordic_angle_event,
	cordic_asine_event,
	cordic_atan_event,
	cordic_cexp_event,
	cordic_sine_event,
	cordic_last_event
} cordic_event;

/* Error_type */
typedef enum
{
	cordic_acosine_errno,
	cordic_cosine_errno,
	cordic_sincos_errno,
	cordic_angle_errno,
	cordic_asine_errno,
	cordic_atan_errno,
	cordic_cexp_errno,
	cordic_sine_errno,
	cordic_last_errno
} eSystemErrno;

enum  ErrNo { CORDIC_ACOSINE, CORDIC_COSINE, CORDIC_SINCOS, CORDIC_ANGLE, CORDIC_ASINE, CORDIC_ATAN, CORDIC_CEXP, CORDIC_SINE, CORDIC_LAST };
/* typedef of 2dim array */
typedef cordicstate(* const tb_cdc_eventhandler[cordic_last_state][cordic_last_event])(void);
/* typedef of function pointer */
typedef cordicstate(*pfEventHandler)(void);
/*  API -cmocka infra to control indexing */
int32_t read_fsm_handler(int32_t idx);
/* API - Errno */
int32_t fsmreaderr(int32_t idx);

/* Function Declarations */
extern void acos_fixpt(const int32_t cdc_acos_th[21], int32_t th_rad_fxp[21]);
extern void asin_fixpt(const int32_t a[21], int32_t th_rad_fxp[21]);
int32_t iScalarCordicACos(int32_t cosValue, int16_t nIters, const int32_t LUT[30]);
int32_t iScalarCordicAsin(int32_t sinValue, int16_t nIters, const int32_t LUT[30]);
extern void init_acos_fixpt(int32_t a[21]);
extern void init_asin_fixpt(int32_t a[21]);
extern void init_cosine_fixpt(int32_t thRadFxp[2047]);
extern void cosine_fixpt(const int32_t thRadFxp[2047], int32_t cdcCoSinTh[2047]);
extern void init_sine_fixpt(int32_t thRadFxp[2047]);
extern void sine_fixpt(const int32_t thRadFxp[2047], int32_t cdcCoSinTh[2047]);
extern void init_atan2_fixpt(int32_t thetax[359], int32_t thetay[359], uint16_t z[359]);
extern void atan2_fixpt(const int32_t thetaY[359], const int32_t thetaX[359],
			int32_t cdcatan2Th[359]);
extern void init_cmpx_exp_fixpt(int32_t thRadFxp[2047]);
extern void cmpx_exp_fixpt(const int32_t thRadFxp[2047], cint32_T cdcCoSinTh[2047]);
extern void init_sincos_fixpt(int32_t theta[2047]);
extern void sincos_fixpt(const int32_t theta[2047], int32_t y[2047], int32_t x[2047]);
extern void init_angle_fixpt(cint32_T dblRandomVals[20]);
extern void angle_fixpt(const cint32_T dblRandomVals[20], int32_t theta_dbl_cdc[20]);



#endif