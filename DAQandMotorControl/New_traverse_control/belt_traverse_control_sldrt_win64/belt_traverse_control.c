/*
 * belt_traverse_control.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "belt_traverse_control".
 *
 * Model version              : 1.29
 * Simulink Coder version : 9.7 (R2022a) 13-Nov-2021
 * C source code generated on : Fri Feb 17 11:53:04 2023
 *
 * Target selection: sldrt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "belt_traverse_control.h"
#include "rtwtypes.h"
#include "rt_nonfinite.h"
#include "belt_traverse_control_dt.h"
#define belt_traverse_control_period   (0.001)

/* options for Simulink Desktop Real-Time board 0 */
static double SLDRTBoardOptions0[] = {
  0.0,
  3.0,
  61440.0,
  3.0,
  3.0,
  3.0,
  3.0,
};

/* list of Simulink Desktop Real-Time timers */
const int SLDRTTimerCount = 1;
const double SLDRTTimers[2] = {
  0.001, 0.0,
};

/* list of Simulink Desktop Real-Time boards */
const int SLDRTBoardCount = 1;
SLDRTBOARD SLDRTBoards[1] = {
  { "National_Instruments/PCIe-6323", 4294967295U, 7, SLDRTBoardOptions0 },
};

/* Block states (default storage) */
DW_belt_traverse_control_T belt_traverse_control_DW;

/* Real-time model */
static RT_MODEL_belt_traverse_control_T belt_traverse_control_M_;
RT_MODEL_belt_traverse_control_T *const belt_traverse_control_M =
  &belt_traverse_control_M_;

/* Model output function */
void belt_traverse_control_output(void)
{
  real_T rtb_Saturation1;
  real_T rtb_Saturation2;
  real_T rtb_Saturation3;
  real_T rtb_Saturation4;
  real_T rtb_theta_cmd_ai1;
  real32_T rtb_theta_lock_in_p03;
  boolean_T rtb_LogicalOperator;
  boolean_T rtb_LogicalOperator1;

  /* S-Function (sldrtdi): '<Root>/theta_lock_in_p0.3' */
  /* S-Function Block: <Root>/theta_lock_in_p0.3 */
  {
    double inval[1];
    double* invalp = inval;
    RTBIO_DriverIO(0, DIGITALINPUT, IOREAD, 1,
                   &belt_traverse_control_P.theta_lock_in_p03_Channels, inval,
                   &belt_traverse_control_P.theta_lock_in_p03_BitMode);
    rtb_theta_lock_in_p03 = (real32_T) *invalp++;
  }

  /* S-Function (sldrtdi): '<Root>/y_lock_in_p0.2' */
  /* S-Function Block: <Root>/y_lock_in_p0.2 */
  {
    double inval[1];
    double* invalp = inval;
    RTBIO_DriverIO(0, DIGITALINPUT, IOREAD, 1,
                   &belt_traverse_control_P.y_lock_in_p02_Channels, inval,
                   &belt_traverse_control_P.y_lock_in_p02_BitMode);
    rtb_LogicalOperator1 = (boolean_T) *invalp++;
  }

  /* S-Function (sldrtai): '<Root>/y_cmd_ai0' */
  /* S-Function Block: <Root>/y_cmd_ai0 */
  {
    ANALOGIOPARM parm;
    parm.mode = (RANGEMODE) belt_traverse_control_P.y_cmd_ai0_RangeMode;
    parm.rangeidx = belt_traverse_control_P.y_cmd_ai0_VoltRange;
    RTBIO_DriverIO(0, ANALOGINPUT, IOREAD, 1,
                   &belt_traverse_control_P.y_cmd_ai0_Channels, &rtb_Saturation4,
                   &parm);
  }

  /* S-Function (sldrtai): '<Root>/theta_cmd_ai1' */
  /* S-Function Block: <Root>/theta_cmd_ai1 */
  {
    ANALOGIOPARM parm;
    parm.mode = (RANGEMODE) belt_traverse_control_P.theta_cmd_ai1_RangeMode;
    parm.rangeidx = belt_traverse_control_P.theta_cmd_ai1_VoltRange;
    RTBIO_DriverIO(0, ANALOGINPUT, IOREAD, 1,
                   &belt_traverse_control_P.theta_cmd_ai1_Channels,
                   &rtb_theta_cmd_ai1, &parm);
  }

  /* Logic: '<Root>/Logical Operator' incorporates:
   *  Constant: '<Root>/theta_lock_manual'
   */
  rtb_LogicalOperator = (belt_traverse_control_P.theta_lock_manual_Value ||
    (rtb_theta_lock_in_p03 != 0.0F));

  /* Logic: '<Root>/Logical Operator1' incorporates:
   *  Constant: '<Root>/y_lock_manual'
   */
  rtb_LogicalOperator1 = (belt_traverse_control_P.y_lock_manual_Value ||
    rtb_LogicalOperator1);

  /* Saturate: '<Root>/Saturation3' */
  if (rtb_Saturation4 > belt_traverse_control_P.Saturation3_UpperSat) {
    rtb_Saturation3 = belt_traverse_control_P.Saturation3_UpperSat;
  } else if (rtb_Saturation4 < belt_traverse_control_P.Saturation3_LowerSat) {
    rtb_Saturation3 = belt_traverse_control_P.Saturation3_LowerSat;
  } else {
    rtb_Saturation3 = rtb_Saturation4;
  }

  /* End of Saturate: '<Root>/Saturation3' */

  /* RateLimiter: '<Root>/Rate Limiter3' incorporates:
   *  Constant: '<Root>/y_offset [m]'
   */
  rtb_Saturation2 = belt_traverse_control_P.y_offsetm_Value -
    belt_traverse_control_DW.PrevY;
  if (rtb_Saturation2 > belt_traverse_control_P.RateLimiter3_RisingLim *
      belt_traverse_control_period) {
    rtb_Saturation1 = belt_traverse_control_P.RateLimiter3_RisingLim *
      belt_traverse_control_period + belt_traverse_control_DW.PrevY;
  } else if (rtb_Saturation2 < belt_traverse_control_P.RateLimiter3_FallingLim *
             belt_traverse_control_period) {
    rtb_Saturation1 = belt_traverse_control_P.RateLimiter3_FallingLim *
      belt_traverse_control_period + belt_traverse_control_DW.PrevY;
  } else {
    rtb_Saturation1 = belt_traverse_control_P.y_offsetm_Value;
  }

  belt_traverse_control_DW.PrevY = rtb_Saturation1;

  /* End of RateLimiter: '<Root>/Rate Limiter3' */

  /* Bias: '<Root>/y_calib_bias [Hz]' incorporates:
   *  Gain: '<Root>/y_calib [Hz//m]'
   *  Gain: '<Root>/y_cmd_V2m'
   *  Sum: '<Root>/Sum2'
   */
  rtb_Saturation1 = (belt_traverse_control_P.y_cmd_V2m_Gain * rtb_Saturation4 +
                     rtb_Saturation1) * belt_traverse_control_P.y_calibHzm_Gain
    + belt_traverse_control_P.y_calib_biasHz_Bias;

  /* Saturate: '<Root>/Saturation1' */
  if (rtb_Saturation1 > belt_traverse_control_P.Saturation1_UpperSat) {
    rtb_Saturation1 = belt_traverse_control_P.Saturation1_UpperSat;
  } else if (rtb_Saturation1 < belt_traverse_control_P.Saturation1_LowerSat) {
    rtb_Saturation1 = belt_traverse_control_P.Saturation1_LowerSat;
  }

  /* End of Saturate: '<Root>/Saturation1' */

  /* Saturate: '<Root>/Saturation4' */
  if (rtb_theta_cmd_ai1 > belt_traverse_control_P.Saturation4_UpperSat) {
    rtb_Saturation4 = belt_traverse_control_P.Saturation4_UpperSat;
  } else if (rtb_theta_cmd_ai1 < belt_traverse_control_P.Saturation4_LowerSat) {
    rtb_Saturation4 = belt_traverse_control_P.Saturation4_LowerSat;
  } else {
    rtb_Saturation4 = rtb_theta_cmd_ai1;
  }

  /* End of Saturate: '<Root>/Saturation4' */

  /* RateLimiter: '<Root>/Rate Limiter6' incorporates:
   *  Constant: '<Root>/theta_offset [deg]'
   */
  rtb_Saturation2 = belt_traverse_control_P.theta_offsetdeg_Value -
    belt_traverse_control_DW.PrevY_d;
  if (rtb_Saturation2 > belt_traverse_control_P.RateLimiter6_RisingLim *
      belt_traverse_control_period) {
    rtb_Saturation2 = belt_traverse_control_P.RateLimiter6_RisingLim *
      belt_traverse_control_period + belt_traverse_control_DW.PrevY_d;
  } else if (rtb_Saturation2 < belt_traverse_control_P.RateLimiter6_FallingLim *
             belt_traverse_control_period) {
    rtb_Saturation2 = belt_traverse_control_P.RateLimiter6_FallingLim *
      belt_traverse_control_period + belt_traverse_control_DW.PrevY_d;
  } else {
    rtb_Saturation2 = belt_traverse_control_P.theta_offsetdeg_Value;
  }

  belt_traverse_control_DW.PrevY_d = rtb_Saturation2;

  /* End of RateLimiter: '<Root>/Rate Limiter6' */

  /* Bias: '<Root>/theta_calib_bias [Hz]' incorporates:
   *  Gain: '<Root>/theta_calib [Hz//deg]'
   *  Gain: '<Root>/theta_cmd_V2deg'
   *  Sum: '<Root>/Sum3'
   */
  rtb_Saturation2 = (belt_traverse_control_P.theta_cmd_V2deg_Gain *
                     rtb_theta_cmd_ai1 + rtb_Saturation2) *
    belt_traverse_control_P.theta_calibHzdeg_Gain +
    belt_traverse_control_P.theta_calib_biasHz_Bias;

  /* Saturate: '<Root>/Saturation2' */
  if (rtb_Saturation2 > belt_traverse_control_P.Saturation2_UpperSat) {
    rtb_Saturation2 = belt_traverse_control_P.Saturation2_UpperSat;
  } else if (rtb_Saturation2 < belt_traverse_control_P.Saturation2_LowerSat) {
    rtb_Saturation2 = belt_traverse_control_P.Saturation2_LowerSat;
  }

  /* End of Saturate: '<Root>/Saturation2' */

  /* S-Function (sldrtdo): '<Root>/theta_lock_out_signal_p0.1' */
  /* S-Function Block: <Root>/theta_lock_out_signal_p0.1 */
  {
    double doval[1];
    double* dovalp = doval;
    *dovalp++ = (double) rtb_LogicalOperator;
    RTBIO_DriverIO(0, DIGITALOUTPUT, IOWRITE, 1,
                   &belt_traverse_control_P.theta_lock_out_signal_p01_Channels,
                   doval,
                   &belt_traverse_control_P.theta_lock_out_signal_p01_BitMode);
  }

  /* S-Function (sldrtdo): '<Root>/y_lock_out_signal_p0.0' */
  /* S-Function Block: <Root>/y_lock_out_signal_p0.0 */
  {
    double doval[1];
    double* dovalp = doval;
    *dovalp++ = (double) rtb_LogicalOperator1;
    RTBIO_DriverIO(0, DIGITALOUTPUT, IOWRITE, 1,
                   &belt_traverse_control_P.y_lock_out_signal_p00_Channels,
                   doval, &belt_traverse_control_P.y_lock_out_signal_p00_BitMode);
  }

  /* S-Function (sldrtao): '<Root>/y_cmd_out_signal_ao0' */
  /* S-Function Block: <Root>/y_cmd_out_signal_ao0 */
  {
    {
      ANALOGIOPARM parm;
      parm.mode = (RANGEMODE)
        belt_traverse_control_P.y_cmd_out_signal_ao0_RangeMode;
      parm.rangeidx = belt_traverse_control_P.y_cmd_out_signal_ao0_VoltRange;
      RTBIO_DriverIO(0, ANALOGOUTPUT, IOWRITE, 1,
                     &belt_traverse_control_P.y_cmd_out_signal_ao0_Channels,
                     ((real_T*) (&rtb_Saturation3)), &parm);
    }
  }

  /* S-Function (sldrtqo): '<Root>/y_control_signal_ctr0_out_PFI12' */
  /* S-Function Block: <Root>/y_control_signal_ctr0_out_PFI12 */
  {
    {
      RTBIO_DriverIO(0, FREQUENCYOUTPUT, IOWRITE, 1,
                     &belt_traverse_control_P.y_control_signal_ctr0_out_PFI12_Channels,
                     ((real_T*) (&rtb_Saturation1)),
                     &belt_traverse_control_P.y_control_signal_ctr0_out_PFI12_Duty);
    }
  }

  /* S-Function (sldrtao): '<Root>/theta_cmd_out_signal_ao1' */
  /* S-Function Block: <Root>/theta_cmd_out_signal_ao1 */
  {
    {
      ANALOGIOPARM parm;
      parm.mode = (RANGEMODE)
        belt_traverse_control_P.theta_cmd_out_signal_ao1_RangeMode;
      parm.rangeidx = belt_traverse_control_P.theta_cmd_out_signal_ao1_VoltRange;
      RTBIO_DriverIO(0, ANALOGOUTPUT, IOWRITE, 1,
                     &belt_traverse_control_P.theta_cmd_out_signal_ao1_Channels,
                     ((real_T*) (&rtb_Saturation4)), &parm);
    }
  }

  /* S-Function (sldrtqo): '<Root>/theta_control_signal_ctr1_PFI13' */
  /* S-Function Block: <Root>/theta_control_signal_ctr1_PFI13 */
  {
    {
      RTBIO_DriverIO(0, FREQUENCYOUTPUT, IOWRITE, 1,
                     &belt_traverse_control_P.theta_control_signal_ctr1_PFI13_Channels,
                     ((real_T*) (&rtb_Saturation2)),
                     &belt_traverse_control_P.theta_control_signal_ctr1_PFI13_Duty);
    }
  }
}

/* Model update function */
void belt_traverse_control_update(void)
{
  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++belt_traverse_control_M->Timing.clockTick0)) {
    ++belt_traverse_control_M->Timing.clockTickH0;
  }

  belt_traverse_control_M->Timing.t[0] =
    belt_traverse_control_M->Timing.clockTick0 *
    belt_traverse_control_M->Timing.stepSize0 +
    belt_traverse_control_M->Timing.clockTickH0 *
    belt_traverse_control_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
void belt_traverse_control_initialize(void)
{
  /* Start for S-Function (sldrtdo): '<Root>/theta_lock_out_signal_p0.1' */

  /* S-Function Block: <Root>/theta_lock_out_signal_p0.1 */

  /* no initial value should be set */

  /* Start for S-Function (sldrtdo): '<Root>/y_lock_out_signal_p0.0' */

  /* S-Function Block: <Root>/y_lock_out_signal_p0.0 */

  /* no initial value should be set */

  /* Start for S-Function (sldrtao): '<Root>/y_cmd_out_signal_ao0' */

  /* S-Function Block: <Root>/y_cmd_out_signal_ao0 */

  /* no initial value should be set */

  /* Start for S-Function (sldrtqo): '<Root>/y_control_signal_ctr0_out_PFI12' */

  /* S-Function Block: <Root>/y_control_signal_ctr0_out_PFI12 */

  /* no initial value should be set */

  /* Start for S-Function (sldrtao): '<Root>/theta_cmd_out_signal_ao1' */

  /* S-Function Block: <Root>/theta_cmd_out_signal_ao1 */

  /* no initial value should be set */

  /* Start for S-Function (sldrtqo): '<Root>/theta_control_signal_ctr1_PFI13' */

  /* S-Function Block: <Root>/theta_control_signal_ctr1_PFI13 */

  /* no initial value should be set */

  /* InitializeConditions for RateLimiter: '<Root>/Rate Limiter3' */
  belt_traverse_control_DW.PrevY = belt_traverse_control_P.RateLimiter3_IC;

  /* InitializeConditions for RateLimiter: '<Root>/Rate Limiter6' */
  belt_traverse_control_DW.PrevY_d = belt_traverse_control_P.RateLimiter6_IC;
}

/* Model terminate function */
void belt_traverse_control_terminate(void)
{
  /* Terminate for S-Function (sldrtdo): '<Root>/theta_lock_out_signal_p0.1' */

  /* S-Function Block: <Root>/theta_lock_out_signal_p0.1 */

  /* no final value should be set */

  /* Terminate for S-Function (sldrtdo): '<Root>/y_lock_out_signal_p0.0' */

  /* S-Function Block: <Root>/y_lock_out_signal_p0.0 */

  /* no final value should be set */

  /* Terminate for S-Function (sldrtao): '<Root>/y_cmd_out_signal_ao0' */

  /* S-Function Block: <Root>/y_cmd_out_signal_ao0 */

  /* no final value should be set */

  /* Terminate for S-Function (sldrtqo): '<Root>/y_control_signal_ctr0_out_PFI12' */

  /* S-Function Block: <Root>/y_control_signal_ctr0_out_PFI12 */

  /* no final value should be set */

  /* Terminate for S-Function (sldrtao): '<Root>/theta_cmd_out_signal_ao1' */

  /* S-Function Block: <Root>/theta_cmd_out_signal_ao1 */

  /* no final value should be set */

  /* Terminate for S-Function (sldrtqo): '<Root>/theta_control_signal_ctr1_PFI13' */

  /* S-Function Block: <Root>/theta_control_signal_ctr1_PFI13 */

  /* no final value should be set */
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  belt_traverse_control_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  belt_traverse_control_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  belt_traverse_control_initialize();
}

void MdlTerminate(void)
{
  belt_traverse_control_terminate();
}

/* Registration function */
RT_MODEL_belt_traverse_control_T *belt_traverse_control(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)belt_traverse_control_M, 0,
                sizeof(RT_MODEL_belt_traverse_control_T));

  /* Initialize timing info */
  {
    int_T *mdlTsMap = belt_traverse_control_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "belt_traverse_control_M points to
       static memory which is guaranteed to be non-NULL" */
    belt_traverse_control_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    belt_traverse_control_M->Timing.sampleTimes =
      (&belt_traverse_control_M->Timing.sampleTimesArray[0]);
    belt_traverse_control_M->Timing.offsetTimes =
      (&belt_traverse_control_M->Timing.offsetTimesArray[0]);

    /* task periods */
    belt_traverse_control_M->Timing.sampleTimes[0] = (0.001);

    /* task offsets */
    belt_traverse_control_M->Timing.offsetTimes[0] = (0.0);
  }

  rtmSetTPtr(belt_traverse_control_M, &belt_traverse_control_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = belt_traverse_control_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    belt_traverse_control_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(belt_traverse_control_M, -1);
  belt_traverse_control_M->Timing.stepSize0 = 0.001;

  /* External mode info */
  belt_traverse_control_M->Sizes.checksums[0] = (1809890267U);
  belt_traverse_control_M->Sizes.checksums[1] = (754553682U);
  belt_traverse_control_M->Sizes.checksums[2] = (3486107148U);
  belt_traverse_control_M->Sizes.checksums[3] = (3445440987U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    belt_traverse_control_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(belt_traverse_control_M->extModeInfo,
      &belt_traverse_control_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(belt_traverse_control_M->extModeInfo,
                        belt_traverse_control_M->Sizes.checksums);
    rteiSetTPtr(belt_traverse_control_M->extModeInfo, rtmGetTPtr
                (belt_traverse_control_M));
  }

  belt_traverse_control_M->solverInfoPtr = (&belt_traverse_control_M->solverInfo);
  belt_traverse_control_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&belt_traverse_control_M->solverInfo, 0.001);
  rtsiSetSolverMode(&belt_traverse_control_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* parameters */
  belt_traverse_control_M->defaultParam = ((real_T *)&belt_traverse_control_P);

  /* states (dwork) */
  belt_traverse_control_M->dwork = ((void *) &belt_traverse_control_DW);
  (void) memset((void *)&belt_traverse_control_DW, 0,
                sizeof(DW_belt_traverse_control_T));

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    belt_traverse_control_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 22;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  belt_traverse_control_M->Sizes.numContStates = (0);/* Number of continuous states */
  belt_traverse_control_M->Sizes.numY = (0);/* Number of model outputs */
  belt_traverse_control_M->Sizes.numU = (0);/* Number of model inputs */
  belt_traverse_control_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  belt_traverse_control_M->Sizes.numSampTimes = (1);/* Number of sample times */
  belt_traverse_control_M->Sizes.numBlocks = (30);/* Number of blocks */
  belt_traverse_control_M->Sizes.numBlockIO = (0);/* Number of block outputs */
  belt_traverse_control_M->Sizes.numBlockPrms = (68);/* Sum of parameter "widths" */
  return belt_traverse_control_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
