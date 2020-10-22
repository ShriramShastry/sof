// SPDX-License-Identifier: BSD-3-Clause
/*
 * File: impnse.c
 *
 * 
 * Author : Sriram Shastry & Seppo Ingalsuo
 */

#include <sof/audio/buffer.h>
#include <sof/audio/component.h>
#include <sof/audio/format.h>
#include <sof/audio/pipeline.h>
#include <sof/audio/impulse_noise/impnse.h>
#include <sof/common.h>
#include <sof/debug/panic.h>
#include <sof/drivers/ipc.h>
#include <sof/lib/alloc.h>
#include <sof/lib/memory.h>
#include <sof/lib/uuid.h>
#include <sof/list.h>
#include <sof/platform.h>
#include <sof/string.h>
#include <sof/trace/trace.h>
#include <ipc/control.h>
#include <ipc/stream.h>
#include <ipc/topology.h>
#include <user/trace.h>
#include <errno.h>
#include <stddef.h>
#include <stdint.h>

static const struct comp_driver comp_impnse;

/* b809efaf-5681-42b1-9ed6-04bb012dd384 */
DECLARE_SOF_RT_UUID("impnse", impnse_uuid, 0xb809efaf, 0x5681, 0x42b1,
		 0x9e, 0xd6, 0x04, 0xbb, 0x01, 0x2d, 0xd3, 0x84);

DECLARE_TR_CTX(impnse_tr, SOF_UUID(impnse_uuid), LOG_LEVEL_INFO);

/**
 * \brief Sets the impulse noise filter in pass through mode.
 * The frequency response of a DCB filter is:
 * H(z) = (1 - z^-1)/(1-Rz^-1).
 * Setting R to 1 makes the filter act as a passthrough component.
 */
static void impnse_set_passthrough(struct comp_data *cd)
{
	comp_cl_info(&comp_impnse, "impnse_set_passthrough()");
	int i;

	for (i = 0; i < PLATFORM_MAX_CHANNELS; i++)
		cd->ImpnseIsOut.f[i] = ONE_Q2_30;
}

/**
 * \brief Initializes the state of the impulse noise Filter
 */
static void impnse_init_state(struct comp_data *cd)
{
	int i;

	for (i = 0; i < PLATFORM_MAX_CHANNELS; i++) {
		cd->state[i].y_prev = 0;
		cd->state[i].x_prev = 0;
	}
}

/**
 * \brief Creates impulse noise Filter component.
 * \return Pointer to impulse noise Filter component device.
 */
static struct comp_dev *impnse_new(const struct comp_driver *drv,
				    struct sof_ipc_comp *comp)
{
	struct comp_dev *dev;
	struct comp_data *cd;
	struct sof_ipc_comp_process *impnse;
	struct sof_ipc_comp_process *ipc_impnse =
		(struct sof_ipc_comp_process *)comp;
	size_t bs = ipc_impnse->size;
	int ret;

	comp_cl_info(&comp_impnse, "impnse_new()");

	dev = comp_alloc(drv, COMP_SIZE(struct sof_ipc_comp_process));
	if (!dev)
		return NULL;

	impnse = COMP_GET_IPC(dev, sof_ipc_comp_process);
	ret = memcpy_s(impnse, sizeof(*impnse), ipc_impnse,
		       sizeof(struct sof_ipc_comp_process));
	assert(!ret);

	cd = rzalloc(SOF_MEM_ZONE_RUNTIME, 0, SOF_MEM_CAPS_RAM, sizeof(*cd));
	if (!cd) {
		rfree(dev);
		return NULL;
	}

	comp_set_drvdata(dev, cd);

	cd->impnse_func = NULL;
	/**
	 * Copy over the coefficients from the blob to cd->R_coeffs
	 * Set passthrough if the size of the blob is invalid
	 */
	if (bs == sizeof(cd->ImpnseIsOut.f)) {
		ret = memcpy_s(cd->ImpnseIsOut.f, bs, ipc_impnse->data, bs);
		assert(!ret);
	} else {
		if (bs > 0)
			comp_cl_warn(&comp_impnse, "impnse_new(), binary blob size %i, expected %i",
				     bs, sizeof(cd->ImpnseIsOut.f));
		impnse_set_passthrough(cd);
	}

	dev->state = COMP_STATE_READY;
	return dev;
}

/**
 * \brief Frees impulse noise Filter component.
 * \param[in,out] dev impulse noise Filter base component device.
 */
static void impnse_free(struct comp_dev *dev)
{
	struct comp_data *cd = comp_get_drvdata(dev);

	comp_info(dev, "impnse_free()");
	rfree(cd);
	rfree(dev);
}

static int impnse_verify_params(struct comp_dev *dev,
				 struct sof_ipc_stream_params *params)
{
	int ret;

	comp_dbg(dev, "impnse_verify_params()");

	ret = comp_verify_params(dev, 0, params);
	if (ret < 0) {
		comp_err(dev, "impnse_verify_params() error: comp_verify_params() failed.");
		return ret;
	}
	
	return 0;
}

/**
 * \brief Sets impulse noise Filter component audio stream parameters.
 * \param[in,out] dev impulse noise Filter base component device.
 * \return Error code.
 *
 * All done in prepare() since we need to know source and sink component params.
 */
static int impnse_params(struct comp_dev *dev,
			  struct sof_ipc_stream_params *params)
{
	int err;

	comp_dbg(dev, "impnse_params()");

	err = impnse_verify_params(dev, params);
	if (err < 0) {
		comp_err(dev, "impnse_params(): pcm params verification failed");
		return -EINVAL;
	}

	return 0;
}

static int impnse_cmd_get_data(struct comp_dev *dev,
				struct sof_ipc_ctrl_data *cdata,
				size_t max_size)
{
	struct comp_data *cd = comp_get_drvdata(dev);
	size_t resp_size;
	int ret = 0;

	switch (cdata->cmd) {
	case SOF_CTRL_CMD_BINARY:
		comp_info(dev, "impnse_cmd_get_data(), SOF_CTRL_CMD_BINARY");

		/* Copy coefficients back to user space */
		resp_size = sizeof(cd->ImpnseIsOut.f);
		comp_info(dev, "impnse_cmd_get_data(), resp_size %u",
			  resp_size);

		if (resp_size > max_size) {
			comp_err(dev, "response size %i exceeds maximum size %i ",
				 resp_size, max_size);
			ret = -EINVAL;
			break;
		}

		ret = memcpy_s(cdata->data->data, cdata->data->size,
			       cd->ImpnseIsOut.f, resp_size);
		assert(!ret);

		cdata->data->abi = SOF_ABI_VERSION;
		cdata->data->size = resp_size;
		break;
	default:
		comp_err(dev, "impnse_cmd_get_data(), invalid command");
		ret = -EINVAL;
	}

	return ret;
}

static int impnse_cmd_set_data(struct comp_dev *dev,
				struct sof_ipc_ctrl_data *cdata)
{
	struct comp_data *cd = comp_get_drvdata(dev);
	size_t req_size = sizeof(cd->ImpnseIsOut.f);
	int ret = 0;

	switch (cdata->cmd) {
	case SOF_CTRL_CMD_BINARY:
		comp_info(dev, "impnse_cmd_set_data(), SOF_CTRL_CMD_BINARY");

		/* Retrieve the binary controls from the packet */
		ret = memcpy_s(cd->ImpnseIsOut.f, req_size, cdata->data->data,
			       req_size);
		assert(!ret);

		break;
	default:
		comp_err(dev, "impnse_set_data(), invalid command %i",
			 cdata->cmd);
		ret = -EINVAL;
	}

	return ret;
}

/**
 * \brief Handles incoming IPC commands for impulse noise Filter component.
 */
static int impnse_cmd(struct comp_dev *dev, int cmd, void *data,
		       int max_data_size)
{
	struct sof_ipc_ctrl_data *cdata = data;
	int ret = 0;

	comp_info(dev, "impnse_cmd()");

	switch (cmd) {
	case COMP_CMD_SET_DATA:
		ret = impnse_cmd_set_data(dev, cdata);
		break;
	case COMP_CMD_GET_DATA:
		ret = impnse_cmd_get_data(dev, cdata, max_data_size);
		break;
	default:
		comp_err(dev, "impnse_cmd(), invalid command (%i)", cmd);
		ret = -EINVAL;
	}

	return ret;
}

/**
 * \brief Sets impulse noise Filter component state.
 * \param[in,out] dev impulse noise Filter base component device.
 * \param[in] cmd Command type.
 * \return Error code.
 */
static int impnse_trigger(struct comp_dev *dev, int cmd)
{
	comp_info(dev, "impnse_trigger()");

	return comp_set_state(dev, cmd);
}

static void impnse_process(struct comp_dev *dev, struct comp_buffer *source,
			    struct comp_buffer *sink, int frames,
			    uint32_t source_bytes, uint32_t sink_bytes)
{
	struct comp_data *cd = comp_get_drvdata(dev);

	buffer_invalidate(source, source_bytes);

	cd->impnse_func(dev, &source->stream, &sink->stream, frames);

	buffer_writeback(sink, sink_bytes);

	/* calc new free and available */
	comp_update_buffer_consume(source, source_bytes);
	comp_update_buffer_produce(sink, sink_bytes);
}

/**
 * \brief Copies and processes stream data.
 * \param[in,out] dev impulse noise Filter base component device.
 * \return Error code.
 */
static int impnse_copy(struct comp_dev *dev)
{
	struct comp_copy_limits cl;
	struct comp_buffer *sourceb;
	struct comp_buffer *sinkb;

	comp_dbg(dev, "impnse_copy()");

	sourceb = list_first_item(&dev->bsource_list, struct comp_buffer,
				  sink_list);
	sinkb = list_first_item(&dev->bsink_list, struct comp_buffer,
				source_list);

	/* Get source, sink, number of frames etc. to process. */
	comp_get_copy_limits_with_lock(sourceb, sinkb, &cl);

	impnse_process(dev, sourceb, sinkb,
			cl.frames, cl.source_bytes, cl.sink_bytes);

	return 0;
}

/**
 * \brief Prepares impulse noise Filter component for processing.
 * \param[in,out] dev impulse noise Filter base component device.
 * \return Error code.
 */
static int impnse_prepare(struct comp_dev *dev)
{
	struct comp_data *cd = comp_get_drvdata(dev);
	struct sof_ipc_comp_config *config = dev_comp_config(dev);
	struct comp_buffer *sourceb;
	struct comp_buffer *sinkb;
	struct0_T ImpnseIsOut;
	struct0_T AudioSteam;
	uint32_t sink_period_bytes;
	int ret;

	comp_info(dev, "impnse_prepare()");

	ret = comp_set_state(dev, COMP_TRIGGER_PREPARE);
	if (ret < 0)
		return ret;

	if (ret == COMP_STATUS_STATE_ALREADY_SET)
		return PPL_STATUS_PATH_STOP;

	/* DC Filter component will only ever have one source and sink buffer */
	sourceb = list_first_item(&dev->bsource_list,
				  struct comp_buffer, sink_list);
	sinkb = list_first_item(&dev->bsink_list,
				struct comp_buffer, source_list);

	/* get source data format */
	cd->source_format = sourceb->stream.frame_fmt;

	/* get sink data format and period bytes */
	cd->sink_format = sinkb->stream.frame_fmt;
	sink_period_bytes =
			audio_stream_period_bytes(&sinkb->stream, dev->frames);

	if (sinkb->stream.size < config->periods_sink * sink_period_bytes) {
		comp_err(dev, "impnse_prepare(), sink buffer size %d is insufficient",
			 sinkb->stream.size);
		ret = -ENOMEM;
		goto err;
	}
	
	impnse_init_state(cd);
    ImpnseIsOut = init_struc_fixpt();         // function initialization - this is  required for wrapper unit testing
    AudioSteam = impnse_fixpt(ImpnseIsOut);   // function initialization - this is  required for impulse noise detection and cancellation
	
	cd->impnse_func = impnse_find_func(cd->source_format);
	if (!cd->impnse_func) {
		comp_err(dev, "impnse_prepare(), No processing function matching frames format");
		ret = -EINVAL;
		goto err;
	}

	comp_info(dev, "impnse_prepare(), source_format=%d, sink_format=%d",
		  cd->source_format, cd->sink_format);

	return 0;

err:
	comp_set_state(dev, COMP_TRIGGER_RESET);
	return ret;
}

/**
 * \brief Resets impulse noise Filter component.
 * \param[in,out] dev impulse noise Filter base component device.
 * \return Error code.
 */
static int impnse_reset(struct comp_dev *dev)
{
	comp_info(dev, "impnse_reset()");

	comp_set_state(dev, COMP_TRIGGER_RESET);
	return 0;
}

/** \brief impulse noise Filter component definition. */
static const struct comp_driver comp_impnse = {
	.type = SOF_COMP_IMPNSE,
	.uid  = SOF_RT_UUID(impnse_uuid),
	.tctx = &impnse_tr,
	.ops  = {
		 .create	 = impnse_new,
		 .free	 = impnse_free,
		 .params	 = impnse_params,
		 .cmd		 = impnse_cmd,
		 .trigger = impnse_trigger,
		 .copy	 = impnse_copy,
		 .prepare = impnse_prepare,
		 .reset	 = impnse_reset,
	},
};

static SHARED_DATA struct comp_driver_info comp_impnse_info = {
	.drv = &comp_impnse,
};

static void sys_comp_impnse_init(void)
{
	comp_register(platform_shared_get(&comp_impnse_info,
					  sizeof(comp_impnse_info)));
}

DECLARE_MODULE(sys_comp_impnse_init);
