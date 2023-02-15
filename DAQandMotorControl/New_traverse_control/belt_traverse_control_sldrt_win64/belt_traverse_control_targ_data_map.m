    function targMap = targDataMap(),

    ;%***********************
    ;% Create Parameter Map *
    ;%***********************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 4;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc paramMap
        ;%
        paramMap.nSections           = nTotSects;
        paramMap.sectIdxOffset       = sectIdxOffset;
            paramMap.sections(nTotSects) = dumSection; %prealloc
        paramMap.nTotData            = -1;

        ;%
        ;% Auto data (belt_traverse_control_P)
        ;%
            section.nData     = 22;
            section.data(22)  = dumData; %prealloc

                    ;% belt_traverse_control_P.y_control_signal_ctr0_out_PFI12_Duty
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% belt_traverse_control_P.theta_control_signal_ctr1_PFI13_Duty
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 1;

                    ;% belt_traverse_control_P.y_lock_in_p02_MaxMissedTicks
                    section.data(3).logicalSrcIdx = 16;
                    section.data(3).dtTransOffset = 2;

                    ;% belt_traverse_control_P.theta_lock_in_p03_MaxMissedTicks
                    section.data(4).logicalSrcIdx = 17;
                    section.data(4).dtTransOffset = 3;

                    ;% belt_traverse_control_P.y_cmd_ai0_MaxMissedTicks
                    section.data(5).logicalSrcIdx = 18;
                    section.data(5).dtTransOffset = 4;

                    ;% belt_traverse_control_P.theta_cmd_ai1_MaxMissedTicks
                    section.data(6).logicalSrcIdx = 19;
                    section.data(6).dtTransOffset = 5;

                    ;% belt_traverse_control_P.y_lock_out_signal_p00_MaxMissedTicks
                    section.data(7).logicalSrcIdx = 20;
                    section.data(7).dtTransOffset = 6;

                    ;% belt_traverse_control_P.theta_lock_out_signal_p01_MaxMissedTicks
                    section.data(8).logicalSrcIdx = 21;
                    section.data(8).dtTransOffset = 7;

                    ;% belt_traverse_control_P.y_cmd_out_signal_ao0_MaxMissedTicks
                    section.data(9).logicalSrcIdx = 22;
                    section.data(9).dtTransOffset = 8;

                    ;% belt_traverse_control_P.y_control_signal_ctr0_out_PFI12_MaxMissedTicks
                    section.data(10).logicalSrcIdx = 23;
                    section.data(10).dtTransOffset = 9;

                    ;% belt_traverse_control_P.theta_cmd_out_signal_ao1_MaxMissedTicks
                    section.data(11).logicalSrcIdx = 24;
                    section.data(11).dtTransOffset = 10;

                    ;% belt_traverse_control_P.theta_control_signal_ctr1_PFI13_MaxMissedTicks
                    section.data(12).logicalSrcIdx = 25;
                    section.data(12).dtTransOffset = 11;

                    ;% belt_traverse_control_P.y_lock_in_p02_YieldWhenWaiting
                    section.data(13).logicalSrcIdx = 26;
                    section.data(13).dtTransOffset = 12;

                    ;% belt_traverse_control_P.theta_lock_in_p03_YieldWhenWaiting
                    section.data(14).logicalSrcIdx = 27;
                    section.data(14).dtTransOffset = 13;

                    ;% belt_traverse_control_P.y_cmd_ai0_YieldWhenWaiting
                    section.data(15).logicalSrcIdx = 28;
                    section.data(15).dtTransOffset = 14;

                    ;% belt_traverse_control_P.theta_cmd_ai1_YieldWhenWaiting
                    section.data(16).logicalSrcIdx = 29;
                    section.data(16).dtTransOffset = 15;

                    ;% belt_traverse_control_P.y_lock_out_signal_p00_YieldWhenWaiting
                    section.data(17).logicalSrcIdx = 30;
                    section.data(17).dtTransOffset = 16;

                    ;% belt_traverse_control_P.theta_lock_out_signal_p01_YieldWhenWaiting
                    section.data(18).logicalSrcIdx = 31;
                    section.data(18).dtTransOffset = 17;

                    ;% belt_traverse_control_P.y_cmd_out_signal_ao0_YieldWhenWaiting
                    section.data(19).logicalSrcIdx = 32;
                    section.data(19).dtTransOffset = 18;

                    ;% belt_traverse_control_P.y_control_signal_ctr0_out_PFI12_YieldWhenWaiting
                    section.data(20).logicalSrcIdx = 33;
                    section.data(20).dtTransOffset = 19;

                    ;% belt_traverse_control_P.theta_cmd_out_signal_ao1_YieldWhenWaiting
                    section.data(21).logicalSrcIdx = 34;
                    section.data(21).dtTransOffset = 20;

                    ;% belt_traverse_control_P.theta_control_signal_ctr1_PFI13_YieldWhenWaiting
                    section.data(22).logicalSrcIdx = 35;
                    section.data(22).dtTransOffset = 21;

            nTotData = nTotData + section.nData;
            paramMap.sections(1) = section;
            clear section

            section.nData     = 22;
            section.data(22)  = dumData; %prealloc

                    ;% belt_traverse_control_P.y_lock_in_p02_BitMode
                    section.data(1).logicalSrcIdx = 36;
                    section.data(1).dtTransOffset = 0;

                    ;% belt_traverse_control_P.theta_lock_in_p03_BitMode
                    section.data(2).logicalSrcIdx = 37;
                    section.data(2).dtTransOffset = 1;

                    ;% belt_traverse_control_P.y_lock_out_signal_p00_BitMode
                    section.data(3).logicalSrcIdx = 38;
                    section.data(3).dtTransOffset = 2;

                    ;% belt_traverse_control_P.theta_lock_out_signal_p01_BitMode
                    section.data(4).logicalSrcIdx = 39;
                    section.data(4).dtTransOffset = 3;

                    ;% belt_traverse_control_P.y_lock_in_p02_Channels
                    section.data(5).logicalSrcIdx = 40;
                    section.data(5).dtTransOffset = 4;

                    ;% belt_traverse_control_P.theta_lock_in_p03_Channels
                    section.data(6).logicalSrcIdx = 41;
                    section.data(6).dtTransOffset = 5;

                    ;% belt_traverse_control_P.y_cmd_ai0_Channels
                    section.data(7).logicalSrcIdx = 42;
                    section.data(7).dtTransOffset = 6;

                    ;% belt_traverse_control_P.theta_cmd_ai1_Channels
                    section.data(8).logicalSrcIdx = 43;
                    section.data(8).dtTransOffset = 7;

                    ;% belt_traverse_control_P.y_lock_out_signal_p00_Channels
                    section.data(9).logicalSrcIdx = 44;
                    section.data(9).dtTransOffset = 8;

                    ;% belt_traverse_control_P.theta_lock_out_signal_p01_Channels
                    section.data(10).logicalSrcIdx = 45;
                    section.data(10).dtTransOffset = 9;

                    ;% belt_traverse_control_P.y_cmd_out_signal_ao0_Channels
                    section.data(11).logicalSrcIdx = 46;
                    section.data(11).dtTransOffset = 10;

                    ;% belt_traverse_control_P.y_control_signal_ctr0_out_PFI12_Channels
                    section.data(12).logicalSrcIdx = 47;
                    section.data(12).dtTransOffset = 11;

                    ;% belt_traverse_control_P.theta_cmd_out_signal_ao1_Channels
                    section.data(13).logicalSrcIdx = 48;
                    section.data(13).dtTransOffset = 12;

                    ;% belt_traverse_control_P.theta_control_signal_ctr1_PFI13_Channels
                    section.data(14).logicalSrcIdx = 49;
                    section.data(14).dtTransOffset = 13;

                    ;% belt_traverse_control_P.y_cmd_ai0_RangeMode
                    section.data(15).logicalSrcIdx = 50;
                    section.data(15).dtTransOffset = 14;

                    ;% belt_traverse_control_P.theta_cmd_ai1_RangeMode
                    section.data(16).logicalSrcIdx = 51;
                    section.data(16).dtTransOffset = 15;

                    ;% belt_traverse_control_P.y_cmd_out_signal_ao0_RangeMode
                    section.data(17).logicalSrcIdx = 52;
                    section.data(17).dtTransOffset = 16;

                    ;% belt_traverse_control_P.theta_cmd_out_signal_ao1_RangeMode
                    section.data(18).logicalSrcIdx = 53;
                    section.data(18).dtTransOffset = 17;

                    ;% belt_traverse_control_P.y_cmd_ai0_VoltRange
                    section.data(19).logicalSrcIdx = 54;
                    section.data(19).dtTransOffset = 18;

                    ;% belt_traverse_control_P.theta_cmd_ai1_VoltRange
                    section.data(20).logicalSrcIdx = 55;
                    section.data(20).dtTransOffset = 19;

                    ;% belt_traverse_control_P.y_cmd_out_signal_ao0_VoltRange
                    section.data(21).logicalSrcIdx = 56;
                    section.data(21).dtTransOffset = 20;

                    ;% belt_traverse_control_P.theta_cmd_out_signal_ao1_VoltRange
                    section.data(22).logicalSrcIdx = 57;
                    section.data(22).dtTransOffset = 21;

            nTotData = nTotData + section.nData;
            paramMap.sections(2) = section;
            clear section

            section.nData     = 22;
            section.data(22)  = dumData; %prealloc

                    ;% belt_traverse_control_P.Saturation3_UpperSat
                    section.data(1).logicalSrcIdx = 58;
                    section.data(1).dtTransOffset = 0;

                    ;% belt_traverse_control_P.Saturation3_LowerSat
                    section.data(2).logicalSrcIdx = 59;
                    section.data(2).dtTransOffset = 1;

                    ;% belt_traverse_control_P.y_offsetm_Value
                    section.data(3).logicalSrcIdx = 60;
                    section.data(3).dtTransOffset = 2;

                    ;% belt_traverse_control_P.RateLimiter3_RisingLim
                    section.data(4).logicalSrcIdx = 61;
                    section.data(4).dtTransOffset = 3;

                    ;% belt_traverse_control_P.RateLimiter3_FallingLim
                    section.data(5).logicalSrcIdx = 62;
                    section.data(5).dtTransOffset = 4;

                    ;% belt_traverse_control_P.RateLimiter3_IC
                    section.data(6).logicalSrcIdx = 63;
                    section.data(6).dtTransOffset = 5;

                    ;% belt_traverse_control_P.y_cmd_V2m_Gain
                    section.data(7).logicalSrcIdx = 64;
                    section.data(7).dtTransOffset = 6;

                    ;% belt_traverse_control_P.y_calibHzm_Gain
                    section.data(8).logicalSrcIdx = 65;
                    section.data(8).dtTransOffset = 7;

                    ;% belt_traverse_control_P.y_calib_biasHz_Bias
                    section.data(9).logicalSrcIdx = 66;
                    section.data(9).dtTransOffset = 8;

                    ;% belt_traverse_control_P.Saturation1_UpperSat
                    section.data(10).logicalSrcIdx = 67;
                    section.data(10).dtTransOffset = 9;

                    ;% belt_traverse_control_P.Saturation1_LowerSat
                    section.data(11).logicalSrcIdx = 68;
                    section.data(11).dtTransOffset = 10;

                    ;% belt_traverse_control_P.Saturation4_UpperSat
                    section.data(12).logicalSrcIdx = 69;
                    section.data(12).dtTransOffset = 11;

                    ;% belt_traverse_control_P.Saturation4_LowerSat
                    section.data(13).logicalSrcIdx = 70;
                    section.data(13).dtTransOffset = 12;

                    ;% belt_traverse_control_P.theta_offsetdeg_Value
                    section.data(14).logicalSrcIdx = 71;
                    section.data(14).dtTransOffset = 13;

                    ;% belt_traverse_control_P.RateLimiter6_RisingLim
                    section.data(15).logicalSrcIdx = 72;
                    section.data(15).dtTransOffset = 14;

                    ;% belt_traverse_control_P.RateLimiter6_FallingLim
                    section.data(16).logicalSrcIdx = 73;
                    section.data(16).dtTransOffset = 15;

                    ;% belt_traverse_control_P.RateLimiter6_IC
                    section.data(17).logicalSrcIdx = 74;
                    section.data(17).dtTransOffset = 16;

                    ;% belt_traverse_control_P.theta_cmd_V2deg_Gain
                    section.data(18).logicalSrcIdx = 75;
                    section.data(18).dtTransOffset = 17;

                    ;% belt_traverse_control_P.theta_calibHzdeg_Gain
                    section.data(19).logicalSrcIdx = 76;
                    section.data(19).dtTransOffset = 18;

                    ;% belt_traverse_control_P.theta_calib_biasHz_Bias
                    section.data(20).logicalSrcIdx = 77;
                    section.data(20).dtTransOffset = 19;

                    ;% belt_traverse_control_P.Saturation2_UpperSat
                    section.data(21).logicalSrcIdx = 78;
                    section.data(21).dtTransOffset = 20;

                    ;% belt_traverse_control_P.Saturation2_LowerSat
                    section.data(22).logicalSrcIdx = 79;
                    section.data(22).dtTransOffset = 21;

            nTotData = nTotData + section.nData;
            paramMap.sections(3) = section;
            clear section

            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% belt_traverse_control_P.y_lock_manual_Value
                    section.data(1).logicalSrcIdx = 80;
                    section.data(1).dtTransOffset = 0;

                    ;% belt_traverse_control_P.theta_lock_manual_Value
                    section.data(2).logicalSrcIdx = 81;
                    section.data(2).dtTransOffset = 1;

            nTotData = nTotData + section.nData;
            paramMap.sections(4) = section;
            clear section


            ;%
            ;% Non-auto Data (parameter)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        paramMap.nTotData = nTotData;



    ;%**************************
    ;% Create Block Output Map *
    ;%**************************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 0;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc sigMap
        ;%
        sigMap.nSections           = nTotSects;
        sigMap.sectIdxOffset       = sectIdxOffset;
        sigMap.nTotData            = -1;

        ;%
        ;% Auto data (belt_traverse_control_B)
        ;%

            ;%
            ;% Non-auto Data (signal)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        sigMap.nTotData = nTotData;



    ;%*******************
    ;% Create DWork Map *
    ;%*******************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 2;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc dworkMap
        ;%
        dworkMap.nSections           = nTotSects;
        dworkMap.sectIdxOffset       = sectIdxOffset;
            dworkMap.sections(nTotSects) = dumSection; %prealloc
        dworkMap.nTotData            = -1;

        ;%
        ;% Auto data (belt_traverse_control_DW)
        ;%
            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% belt_traverse_control_DW.PrevY
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% belt_traverse_control_DW.PrevY_d
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 1;

            nTotData = nTotData + section.nData;
            dworkMap.sections(1) = section;
            clear section

            section.nData     = 10;
            section.data(10)  = dumData; %prealloc

                    ;% belt_traverse_control_DW.y_lock_in_p02_PWORK
                    section.data(1).logicalSrcIdx = 2;
                    section.data(1).dtTransOffset = 0;

                    ;% belt_traverse_control_DW.theta_lock_in_p03_PWORK
                    section.data(2).logicalSrcIdx = 3;
                    section.data(2).dtTransOffset = 1;

                    ;% belt_traverse_control_DW.y_cmd_ai0_PWORK
                    section.data(3).logicalSrcIdx = 4;
                    section.data(3).dtTransOffset = 2;

                    ;% belt_traverse_control_DW.theta_cmd_ai1_PWORK
                    section.data(4).logicalSrcIdx = 5;
                    section.data(4).dtTransOffset = 3;

                    ;% belt_traverse_control_DW.y_lock_out_signal_p00_PWORK
                    section.data(5).logicalSrcIdx = 6;
                    section.data(5).dtTransOffset = 4;

                    ;% belt_traverse_control_DW.theta_lock_out_signal_p01_PWORK
                    section.data(6).logicalSrcIdx = 7;
                    section.data(6).dtTransOffset = 5;

                    ;% belt_traverse_control_DW.y_cmd_out_signal_ao0_PWORK
                    section.data(7).logicalSrcIdx = 8;
                    section.data(7).dtTransOffset = 6;

                    ;% belt_traverse_control_DW.y_control_signal_ctr0_out_PFI12_PWORK
                    section.data(8).logicalSrcIdx = 9;
                    section.data(8).dtTransOffset = 7;

                    ;% belt_traverse_control_DW.theta_cmd_out_signal_ao1_PWORK
                    section.data(9).logicalSrcIdx = 10;
                    section.data(9).dtTransOffset = 8;

                    ;% belt_traverse_control_DW.theta_control_signal_ctr1_PFI13_PWORK
                    section.data(10).logicalSrcIdx = 11;
                    section.data(10).dtTransOffset = 9;

            nTotData = nTotData + section.nData;
            dworkMap.sections(2) = section;
            clear section


            ;%
            ;% Non-auto Data (dwork)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        dworkMap.nTotData = nTotData;



    ;%
    ;% Add individual maps to base struct.
    ;%

    targMap.paramMap  = paramMap;
    targMap.signalMap = sigMap;
    targMap.dworkMap  = dworkMap;

    ;%
    ;% Add checksums to base struct.
    ;%


    targMap.checksum0 = 2861307693;
    targMap.checksum1 = 2657055331;
    targMap.checksum2 = 1646048555;
    targMap.checksum3 = 3412688750;

