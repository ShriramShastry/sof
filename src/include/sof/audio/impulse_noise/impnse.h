 /*
 * SPDX - License - Identifier: BSD - 3 - Clause
  * File: impnse.h
 *
 * Author : Sriram Shastry & Ingalsuo, Seppo
 * Coyright(c) 2017 Intel Corporation.All rights reserved.
 */

#ifndef IMPNSE_H
#define IMPNSE_H

/* Include Files */
#include <stddef.h>
#include <sof/audio/impulse_noise/stdmlib.h>

#include <sof/audio/impulse_noise/impnse_types.h>
#include <stdint.h>
#include <sof/platform.h>
#include <ipc/stream.h>
#include <sof/audio/impulse_noise/rtwtypes.h>

/* Function Declarations */
extern struct0_T impnse_fixpt(const struct0_T ImpnseIsOut);
extern void impnse_initialize(void);
extern void impnse_terminate(void);
extern struct0_T init_struc_fixpt(void);


struct audio_stream;
struct comp_dev;

struct impnse_state {
	int32_t x_prev; /**< state variable referring to x[n-1] */
	int32_t y_prev; /**< state variable referring to y[n-1] */
};


struct impnse_T {
    unsigned char NumChan;
    int f[202];
    uint64m_T minHeight[2];
    uint64m_T repValue[2];
    int64m_T fFilt[202];
    bool isOut[202];
};


/**
 * \brief Type definition for the processing function for the
 * DC Blocking Filter.
 */
typedef void (*impnse_func)(const struct comp_dev *dev,
			     const struct audio_stream *source,
			     const struct audio_stream *sink,
			     uint32_t frames);

/* DC Blocking Filter component private data */
struct comp_data {
	/**< filters state */
	struct impnse_state state[PLATFORM_MAX_CHANNELS];

	/** coefficients for the processing function */
	int32_t R_coeffs[PLATFORM_MAX_CHANNELS];
    struct impnse_T ImpnseIsOut;
	enum sof_ipc_frame source_format;
	enum sof_ipc_frame sink_format;
	impnse_func impnse_func; /**< processing function */
};

/** \brief DC Blocking Filter processing functions map item. */
struct impnse_func_map {
	enum sof_ipc_frame src_fmt; /**< source frame format */
	impnse_func func; /**< processing function */
};

/** \brief Map of formats with dedicated processing functions. */
extern const struct impnse_func_map impnse_fnmap[];

/** \brief Number of processing functions. */
extern const size_t impnse_fncount;

/**
 * \brief Retrieves a DC Blocking processing function matching
 *	  the source buffer's frame format.
 * \param src_fmt the frames' format of the source buffer
 */
static inline impnse_func impnse_find_func(enum sof_ipc_frame src_fmt)
{
	int i;

	/* Find suitable processing function from map */
	for (i = 0; i < impnse_fncount; i++) {
		if (src_fmt == impnse_fnmap[i].src_fmt)
			return impnse_fnmap[i].func;
	}

	return NULL;
}


#endif

/*
 * File trailer for impnse.h
 *
 * [EOF]
 */
