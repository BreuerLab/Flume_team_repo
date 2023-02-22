/*
 * belt_traverse_control_data.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "belt_traverse_control".
 *
 * Model version              : 1.30
 * Simulink Coder version : 9.7 (R2022a) 13-Nov-2021
 * C source code generated on : Tue Feb 21 19:17:28 2023
 *
 * Target selection: sldrt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "belt_traverse_control.h"

/* Block parameters (default storage) */
P_belt_traverse_control_T belt_traverse_control_P = {
  /* Mask Parameter: y_control_signal_ctr0_out_PFI12_Duty
   * Referenced by: '<Root>/y_control_signal_ctr0_out_PFI12'
   */
  0.5,

  /* Mask Parameter: theta_control_signal_ctr1_PFI13_Duty
   * Referenced by: '<Root>/theta_control_signal_ctr1_PFI13'
   */
  0.5,

  /* Mask Parameter: theta_lock_in_p03_MaxMissedTicks
   * Referenced by: '<Root>/theta_lock_in_p0.3'
   */
  10.0,

  /* Mask Parameter: y_lock_in_p02_MaxMissedTicks
   * Referenced by: '<Root>/y_lock_in_p0.2'
   */
  10.0,

  /* Mask Parameter: y_cmd_ai0_MaxMissedTicks
   * Referenced by: '<Root>/y_cmd_ai0'
   */
  10.0,

  /* Mask Parameter: theta_cmd_ai1_MaxMissedTicks
   * Referenced by: '<Root>/theta_cmd_ai1'
   */
  10.0,

  /* Mask Parameter: theta_lock_out_signal_p01_MaxMissedTicks
   * Referenced by: '<Root>/theta_lock_out_signal_p0.1'
   */
  10.0,

  /* Mask Parameter: y_lock_out_signal_p00_MaxMissedTicks
   * Referenced by: '<Root>/y_lock_out_signal_p0.0'
   */
  10.0,

  /* Mask Parameter: y_cmd_out_signal_ao0_MaxMissedTicks
   * Referenced by: '<Root>/y_cmd_out_signal_ao0'
   */
  10.0,

  /* Mask Parameter: y_control_signal_ctr0_out_PFI12_MaxMissedTicks
   * Referenced by: '<Root>/y_control_signal_ctr0_out_PFI12'
   */
  10.0,

  /* Mask Parameter: theta_cmd_out_signal_ao1_MaxMissedTicks
   * Referenced by: '<Root>/theta_cmd_out_signal_ao1'
   */
  10.0,

  /* Mask Parameter: theta_control_signal_ctr1_PFI13_MaxMissedTicks
   * Referenced by: '<Root>/theta_control_signal_ctr1_PFI13'
   */
  10.0,

  /* Mask Parameter: theta_lock_in_p03_YieldWhenWaiting
   * Referenced by: '<Root>/theta_lock_in_p0.3'
   */
  1.0,

  /* Mask Parameter: y_lock_in_p02_YieldWhenWaiting
   * Referenced by: '<Root>/y_lock_in_p0.2'
   */
  1.0,

  /* Mask Parameter: y_cmd_ai0_YieldWhenWaiting
   * Referenced by: '<Root>/y_cmd_ai0'
   */
  0.0,

  /* Mask Parameter: theta_cmd_ai1_YieldWhenWaiting
   * Referenced by: '<Root>/theta_cmd_ai1'
   */
  0.0,

  /* Mask Parameter: theta_lock_out_signal_p01_YieldWhenWaiting
   * Referenced by: '<Root>/theta_lock_out_signal_p0.1'
   */
  0.0,

  /* Mask Parameter: y_lock_out_signal_p00_YieldWhenWaiting
   * Referenced by: '<Root>/y_lock_out_signal_p0.0'
   */
  0.0,

  /* Mask Parameter: y_cmd_out_signal_ao0_YieldWhenWaiting
   * Referenced by: '<Root>/y_cmd_out_signal_ao0'
   */
  0.0,

  /* Mask Parameter: y_control_signal_ctr0_out_PFI12_YieldWhenWaiting
   * Referenced by: '<Root>/y_control_signal_ctr0_out_PFI12'
   */
  0.0,

  /* Mask Parameter: theta_cmd_out_signal_ao1_YieldWhenWaiting
   * Referenced by: '<Root>/theta_cmd_out_signal_ao1'
   */
  0.0,

  /* Mask Parameter: theta_control_signal_ctr1_PFI13_YieldWhenWaiting
   * Referenced by: '<Root>/theta_control_signal_ctr1_PFI13'
   */
  0.0,

  /* Mask Parameter: theta_lock_in_p03_BitMode
   * Referenced by: '<Root>/theta_lock_in_p0.3'
   */
  0,

  /* Mask Parameter: y_lock_in_p02_BitMode
   * Referenced by: '<Root>/y_lock_in_p0.2'
   */
  0,

  /* Mask Parameter: theta_lock_out_signal_p01_BitMode
   * Referenced by: '<Root>/theta_lock_out_signal_p0.1'
   */
  0,

  /* Mask Parameter: y_lock_out_signal_p00_BitMode
   * Referenced by: '<Root>/y_lock_out_signal_p0.0'
   */
  0,

  /* Mask Parameter: theta_lock_in_p03_Channels
   * Referenced by: '<Root>/theta_lock_in_p0.3'
   */
  3,

  /* Mask Parameter: y_lock_in_p02_Channels
   * Referenced by: '<Root>/y_lock_in_p0.2'
   */
  2,

  /* Mask Parameter: y_cmd_ai0_Channels
   * Referenced by: '<Root>/y_cmd_ai0'
   */
  0,

  /* Mask Parameter: theta_cmd_ai1_Channels
   * Referenced by: '<Root>/theta_cmd_ai1'
   */
  1,

  /* Mask Parameter: theta_lock_out_signal_p01_Channels
   * Referenced by: '<Root>/theta_lock_out_signal_p0.1'
   */
  1,

  /* Mask Parameter: y_lock_out_signal_p00_Channels
   * Referenced by: '<Root>/y_lock_out_signal_p0.0'
   */
  0,

  /* Mask Parameter: y_cmd_out_signal_ao0_Channels
   * Referenced by: '<Root>/y_cmd_out_signal_ao0'
   */
  0,

  /* Mask Parameter: y_control_signal_ctr0_out_PFI12_Channels
   * Referenced by: '<Root>/y_control_signal_ctr0_out_PFI12'
   */
  0,

  /* Mask Parameter: theta_cmd_out_signal_ao1_Channels
   * Referenced by: '<Root>/theta_cmd_out_signal_ao1'
   */
  1,

  /* Mask Parameter: theta_control_signal_ctr1_PFI13_Channels
   * Referenced by: '<Root>/theta_control_signal_ctr1_PFI13'
   */
  1,

  /* Mask Parameter: y_cmd_ai0_RangeMode
   * Referenced by: '<Root>/y_cmd_ai0'
   */
  0,

  /* Mask Parameter: theta_cmd_ai1_RangeMode
   * Referenced by: '<Root>/theta_cmd_ai1'
   */
  0,

  /* Mask Parameter: y_cmd_out_signal_ao0_RangeMode
   * Referenced by: '<Root>/y_cmd_out_signal_ao0'
   */
  0,

  /* Mask Parameter: theta_cmd_out_signal_ao1_RangeMode
   * Referenced by: '<Root>/theta_cmd_out_signal_ao1'
   */
  0,

  /* Mask Parameter: y_cmd_ai0_VoltRange
   * Referenced by: '<Root>/y_cmd_ai0'
   */
  0,

  /* Mask Parameter: theta_cmd_ai1_VoltRange
   * Referenced by: '<Root>/theta_cmd_ai1'
   */
  0,

  /* Mask Parameter: y_cmd_out_signal_ao0_VoltRange
   * Referenced by: '<Root>/y_cmd_out_signal_ao0'
   */
  0,

  /* Mask Parameter: theta_cmd_out_signal_ao1_VoltRange
   * Referenced by: '<Root>/theta_cmd_out_signal_ao1'
   */
  0,

  /* Expression: 10
   * Referenced by: '<Root>/Saturation3'
   */
  10.0,

  /* Expression: -10
   * Referenced by: '<Root>/Saturation3'
   */
  -10.0,

  /* Expression: 0
   * Referenced by: '<Root>/y_offset [m]'
   */
  0.0,

  /* Expression: 0.08
   * Referenced by: '<Root>/Rate Limiter3'
   */
  0.08,

  /* Expression: -0.08
   * Referenced by: '<Root>/Rate Limiter3'
   */
  -0.08,

  /* Expression: 0
   * Referenced by: '<Root>/Rate Limiter3'
   */
  0.0,

  /* Expression: 0.05
   * Referenced by: '<Root>/y_cmd_V2m'
   */
  0.05,

  /* Expression: 10000
   * Referenced by: '<Root>/y_calib [Hz//m]'
   */
  10000.0,

  /* Expression: 5000
   * Referenced by: '<Root>/y_calib_bias [Hz]'
   */
  5000.0,

  /* Expression: 10000
   * Referenced by: '<Root>/Saturation1'
   */
  10000.0,

  /* Expression: 5000
   * Referenced by: '<Root>/Saturation1'
   */
  5000.0,

  /* Expression: 10
   * Referenced by: '<Root>/Saturation4'
   */
  10.0,

  /* Expression: -10
   * Referenced by: '<Root>/Saturation4'
   */
  -10.0,

  /* Expression: 0
   * Referenced by: '<Root>/theta_offset [deg]'
   */
  0.0,

  /* Expression: 30
   * Referenced by: '<Root>/Rate Limiter6'
   */
  30.0,

  /* Expression: -30
   * Referenced by: '<Root>/Rate Limiter6'
   */
  -30.0,

  /* Expression: 0
   * Referenced by: '<Root>/Rate Limiter6'
   */
  0.0,

  /* Expression: 36
   * Referenced by: '<Root>/theta_cmd_V2deg'
   */
  36.0,

  /* Expression: 13.8889
   * Referenced by: '<Root>/theta_calib [Hz//deg]'
   */
  13.8889,

  /* Expression: 5000
   * Referenced by: '<Root>/theta_calib_bias [Hz]'
   */
  5000.0,

  /* Expression: 9500
   * Referenced by: '<Root>/Saturation2'
   */
  9500.0,

  /* Expression: 5000
   * Referenced by: '<Root>/Saturation2'
   */
  5000.0,

  /* Computed Parameter: theta_lock_manual_Value
   * Referenced by: '<Root>/theta_lock_manual'
   */
  false,

  /* Computed Parameter: y_lock_manual_Value
   * Referenced by: '<Root>/y_lock_manual'
   */
  false
};
