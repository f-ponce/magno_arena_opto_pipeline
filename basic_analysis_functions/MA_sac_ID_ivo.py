def findSacs(headingVelo,headingVelo_sg, headingVeloThresh = 60.):
    """
    Finds local min-/maxima that exceed a threshold and and sets these saccade epochs to NaNs
    Usage:
    $ import Sac_ID_fp_basic as sid_fp
    $ headingVelo_noSacs = sid_fp.findSacs(flyAngularVeloLowPassFiltered, [headingVeloThresholdDeg])
    IN:
    headingVelo: the Fly Angular velocity vector
    OUT:
    headingVelo_noSacs: fly angular velocity vector with saccades removed (NaNs where there are saccades)
    """
    import numpy as np

    #headingVeloCeil     = 360.
    headingVelo_noSacs  = np.copy(headingVelo)

    def findS_s_e_Idcs(idx, angV): # calculate difference in L-R WSA between nearest local extremes (of opposite sign)
            lastIdx         = angV.shape[0]-1
            # find the startIndex of the saccade
            k               = 0
            onIdx           = None
            chkMore         = True
            refSign         = np.sign(angV[idx])
            while chkMore:
                k += 1
                compSign    = np.sign(angV[idx-k])
                if refSign != compSign: # if not equal
                    onIdx   = idx-k+1
                    chkMore = False
            # find the stopIndex of the saccade
            k               = 0
            offIdx          = None
            chkMore         = True
            while chkMore and idx+k+1 < lastIdx:
                k += 1
                compSign    = np.sign(angV[idx+k])
                if refSign != compSign:
                    offIdx  = idx+k
                    chkMore = False
            if offIdx == None:
                offIdx      = lastIdx
            return onIdx, offIdx

    # In words: find local extremes that exceed headingVelo threshold during flight

    ## Place holders
    headingVeloMxThr,   headingVeloMnThr= np.zeros_like(headingVelo), np.zeros_like(headingVelo)
    sacIdcxAll=np.zeros_like(headingVelo)*np.nan
    sacIdcxl=np.zeros_like(headingVelo)*np.nan
    sacIdcxr=np.zeros_like(headingVelo)*np.nan
    ## Find local extremes in the rate, that exceed the saccade velocity threshold
    reboundIds = 5#15
    for idx in range (reboundIds,len(headingVelo)-reboundIds):

        if headingVelo[idx]>=headingVelo[idx-1] and headingVelo[idx]>=headingVelo[idx+1] and headingVelo[idx]>headingVeloThresh :#and headingVelo[idx]<headingVeloCeil:
            headingVeloMxThr[idx] = headingVelo[idx]
            sacIdcxAll[idx] = headingVelo[idx]
            sacIdcxr[idx] = headingVelo[idx]
        elif headingVelo[idx]<=headingVelo[idx-1] and headingVelo[idx]<=headingVelo[idx+1] and headingVelo[idx]<-headingVeloThresh :#and headingVelo[idx]>-headingVeloCeil:
            headingVeloMnThr[idx] = headingVelo[idx]
            sacIdcxAll[idx] = headingVelo[idx]
            sacIdcxl[idx] = headingVelo[idx]

    ## set saccade epochs to nan
    for idx in range (int(1.4*reboundIds),len(headingVelo)-int(1.4*reboundIds)):
        if abs(headingVeloMxThr[idx])>0: # for every local max
            onIdx, offIdx = findS_s_e_Idcs(idx, headingVelo)
            headingVelo_noSacs [onIdx:offIdx] = np.nan
        if abs(headingVeloMnThr[idx])>0: # for every local minimum
            onIdx, offIdx = findS_s_e_Idcs(idx, headingVelo)
            headingVelo_noSacs [onIdx:offIdx] = np.nan

    return headingVelo_noSacs,sacIdcxAll, sacIdcxl, sacIdcxr
