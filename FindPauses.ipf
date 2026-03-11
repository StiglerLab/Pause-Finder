#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=FindPauses

Menu "Analysis"
	"Find Pauses", FindPausesMenu()
End


Function FindPausesMenu()
	String wylist = WaveList("*", ";", "")
	String wxlist = "__Use x scaling__;"+WaveList("*", ";", "")
	String sigmalist = "__none__;"+WaveList("*", ";", "")
	
	String wnameY, wnameX, wnamesigma
	Prompt wnameY, "Y wave", popup, wylist
	Prompt wnameX, "X wave", popup, wxlist
	Prompt wnameSigma, "Y err (=sigma) wave", popup, sigmalist
	DoPrompt "Find Pauses", wnameY, wnameX, wnameSigma
	
	Wave wy=$wnameY, wx=$wnameX, wsigma=$wnameSigma
	
	if(!WaveExists(wy))
		return NaN
	endif
	
	Variable useXwave = 0
	if(WaveExists(wx))
		if(numpnts(wx) != numpnts(wy))
			abort "X wave and Y wave need to have the same number of points"
		endif
		useXwave = 1
	endif
	if(WaveExists(wSigma))
		if(numpnts(wsigma) != numpnts(wy))
			abort "sigma wave and Y wave need to have the same number of points"
		endif
	endif
	
	Wave W_Pauses = FindPausesMCMC_MT(wy, xwave=wx, sigmawave=wSigma)
	DFREF fld = getwavesdatafolderdfr(W_Pauses)
	Wave/SDFR=fld traj, bestfit, bestfit_colorwave, xwave
	DoWindow FindPausesResults
	if(!v_flag)
		if(useXwave)
			Display/N=FindPausesResults traj, bestfit vs xwave
		else
			Display/N=FindPausesResults traj, bestfit
		endif
		ModifyGraph mode(traj)=2,rgb(traj)=(0,0,0),zColor(bestFit)={:packages:FindPauses:bestFit_colorwave,0,1,RedWhiteGreen,1}
	endif
End



//#define USE_COM //use COM to get end slopes (faster, but does not use varianceWave)

//---- Compute the Bayesian/Schwarz information criterion
//NoiseSigma: Externally provided sigma of noise. If not given, it is estimated (works well for long trajectories only)
//xwave: Optional x-wave with same number of points as traj for non-equal x-spacing.
//varianceWave: expected variance per point (optional). If provided, it precludes noiseSigma. Variance should be tether movement variance + localization variance
Threadsafe Static Function GetSIC(Wave traj, Wave W_pauses [Variable noiseSigma, Wave xwave, Wave varianceWave, Variable sum_ln2piVariances])
	DFREF fld = root:Packages:FindPauses

	Variable useVarianceWave = !paramisdefault(varianceWave) && WaveExists(varianceWave)
	if(useVarianceWave && paramIsDefault(noiseSigma)==0 && numtype(noiseSigma)==0 && noiseSigma!=0)
		print "Cannot use both noiseSigma and variance wave simultaneously!"
		return NaN
	endif
		
	Make/O/N=(numpnts(traj)) fld:fit/Wave=Fit = NaN
	CopyScales/P traj, fit

	Variable useXwave = 0
	if(!paramisdefault(xwave) && waveExists(xwave))
		useXwave = 1
		if(numpnts(xwave) != numpnts(traj))
			print "GetSIC: xwave given but it has different number of points than traj."
			return NaN
		endif
	endif
	
	Duplicate/O W_Pauses, fld:W_PausesUSED/Wave=W_PausesUSED
	FindPauses#SortPauses(W_PausesUSED)
	
	Variable Npauses = dimsize(W_PausesUSED,0)

	//--- Fit pauses
	Variable firstPauseStart = inf, lastPauseEnd = -inf
	Variable i
	Variable v_PausesOK = 0
	for(i=Npauses-1;i>=0;i-=1)
		Variable pauseStart = W_PausesUSED[i][0]
		Variable pauseEnd = pauseStart+W_PausesUSED[i][1]
		
		firstPauseStart = min(firstPauseStart, pauseStart)
		lastPauseEnd = max(lastPauseEnd, pauseEnd)
		
		if(pauseStart < 0 || pauseEnd >= numpnts(traj)) //this shouldn't happen but does (pauseEnd<0)...
			return inf
		endif
		
		Variable wmean = getWeightedMean(pauseStart, pauseEnd)
		fit[pauseStart,pauseEnd] = wmean
		v_pausesOK = v_pausesOK || !numtype(wmean)
	endfor
	
	//--- fit slopes (only first and last are truly "free")
	Variable slope
	if(Npauses==0)
		CurveFit/Q line, traj/X=xwave/NWOK
		Wave w_coef
		if(useXwave)
			fit = w_coef[0] + w_coef[1]*xwave[p]
		else
			fit = w_coef[0] + w_coef[1]*x
		endif
	else
		if(!v_pausesOK)
			return NaN
		endif
		
		
		//---fit before first pause
		Variable x0, y0
		if(firstPauseStart != 0)
			if(useXwave)
				x0 = xwave[firstPauseStart]
			else
				x0 = pnt2x(traj,firstPauseStart)
			endif
			y0 = fit[firstPauseStart]
			
#ifdef USE_COM
			//COM
			Wavestats/Q/M=1/R=[0,firstPauseStart] traj
			Variable COMY = v_avg
			if(!useXwave)
				Duplicate/O/FREE/R=[0,firstPauseStart] traj, wx
				wx = !numtype(traj) ? x : NaN
			else
				Duplicate/O/FREE/R=[0,firstPauseStart] xwave, wx
				wx = !numtype(traj) ? wx[p] : NaN
			endif
			Wavestats/Q/M=1 wx
			Variable COMX = v_avg
			
			slope = (COMY-y0)/(COMX-x0)
#else			
			Wave/SDFR=fld SlopeFromLeftEnd
			slope = SlopeFromLeftEnd[firstPauseStart]
			//slope = getSlope(traj,0,firstPauseStart, xwave=xwave)
#endif

			if(numtype(slope))
				return NaN
			else
				if(!useXwave)
					//Multithread fit[0,firstPauseStart-1] = lineFitOffset({slope, x0, y0},x) //only fill first point: Interpolate2 below takes care of the rest (faster)
					fit[0] = lineFitOffset({slope, x0, y0},pnt2x(fit,0))
				else
					//Multithread fit[0,firstPauseStart-1] = lineFitOffset({slope, x0, y0},xwave[p]) //only fill first point: Interpolate2 below takes care of the rest (faster)
					fit[0] = lineFitOffset({slope, x0, y0},xwave[0])
				endif
			endif
		endif
		

		//---fit after last pause
		if(lastPauseEnd != numpnts(traj)-1)
			if(useXwave)
				x0 = xwave[lastPauseEnd]
			else
				x0 = pnt2x(traj,lastPauseEnd)
			endif
			y0 = fit[lastPauseEnd]
	
#ifdef USE_COM
			//COM
			Wavestats/Q/M=1/R=[lastPauseEnd,] traj
			COMY = v_avg
			if(!useXwave)
				Duplicate/O/FREE/R=[lastPauseEnd,] traj, wx
				wx = !numtype(traj[p+lastPauseEnd]) ? x : NaN
			else
				Duplicate/O/FREE/R=[lastPauseEnd,] xwave, wx
				wx = !numtype(traj[p+lastPauseEnd]) ? wx[p] : NaN
			endif
			Wavestats/Q/M=1 wx
			COMX = v_avg

			slope = (COMY-y0)/(COMX-x0)
#else
			Wave/SDFR=fld SlopeFromRightEnd
			slope = SlopeFromRightEnd[lastPauseEnd]
			//slope = getSlope(traj,lastPauseEnd,numpnts(traj)-1, xwave=xwave)
#endif
			if(numtype(slope))
				return NaN
			else
				if(!useXwave)
					//Multithread fit[lastPauseEnd+1,] = lineFitOffset({slope, x0, y0},x)
					fit[numpnts(fit)-1] = lineFitOffset({slope, x0, y0},x) //only fill last point: Interpolate2 below takes care of the rest (faster)
				else
					//Multithread fit[lastPauseEnd+1,] = lineFitOffset({slope, x0, y0},xwave[p])
					fit[numpnts(fit)-1] = lineFitOffset({slope, x0, y0},xwave[p]) //only fill last point: Interpolate2 below takes care of the rest (faster)
				endif
			endif
		endif
	endif

	if(useXwave)
		Duplicate/O/FREE xwave, xwave_filled //dummy
		Interpolate2/T=1/N=(numpnts(fit))/Y=fit_filled/I=3/X=xwave_filled xwave, fit //fill space between pauses
	else
		Interpolate2/T=1/N=(numpnts(fit))/Y=fit_filled fit //fill space between pauses
	endif
	Duplicate/O fit_filled, fld:fit/Wave=fit
	
	//------ Calculate SIC --------
	//below: numpnts(traj) != n (!!!) 
	if(Npauses>=2)
		Variable Npar = Npauses + 2*(Npauses-1) + 2*(firstPauseStart != 0) + 2*(lastPauseEnd != numpnts(traj)-1) + 1 //1 for mean per pause, 2 DOFs for slopes in between (intercept+slope or times of two adjacent pauses), 2 each for flanking slopes, 1 for global variance
	elseif(Npauses==1)
		Npar = Npauses + 2*(firstPauseStart != 0) + 2*(lastPauseEnd != numpnts(traj)-1) + 1 //1 for pause mean, 1 each for potential flanking lines, 1 global variance
	elseif(Npauses==0)
		Npar = 2 + 1 //two for line fit + 1 for global variance
	else
		print "Error"
		return NaN
	endif

	//Note: The same Npars come out if I compute it like this:
	//starts with slope: y-offset, angle, distance to next fix point =3
	//then distance to fix point of pause +1
	//then angle and distance to next fix point... (if end: distance is not a free parameter)

	if(!useVarianceWave) //---- no variance wave provided (standard). Assumes that every data point has same sigma.
		MatrixOp/O/FREE diff = (traj-fit)

		Wavestats/Q diff
		Variable n = v_npnts
		Variable varMLE = v_sdev^2*(n-1)/n
		Variable/G fld:noiseSigma = sqrt(varMLE)
		Variable/G fld:noiseSigmaEstimated = sqrt(varMLE)
		//Definition: SIC = -2*ln(Likelihood) + Npar * ln(n)
		Variable SIC = (Npar)*ln(n) + n*ln(varMLE) + n*ln(2*pi) + n //SIC
	
		if(!paramIsDefault(noiseSigma) && !numtype(noiseSigma) && noiseSigma!=0) //is noiseSigma provided externally?
			Variable noiseVar = noiseSigma^2 //variance of noise
			Variable/G fld:noiseSigma = noiseSigma
			Npar -= 1 // variance is no longer fit parameter then (?)
			SIC = (Npar)*ln(n) + n*ln(noiseVar) + n*ln(2*pi) + n*(varMLE/noiseVar)
		endif
		
		Variable/G fld:redchisq = n/(n-Npar) * varMLE/noiseVar
	else //---- with variance wave
		MatrixOp/O/FREE scaledResidualsSq = (traj-fit)*(traj-fit)/varianceWave //chance to speed this up by doing cumulative sums of the moments
		Wavestats/Q/M=1 scaledresidualssq
		Variable sum_scaledresidualssq = v_sum
		n = v_npnts
		if(paramisdefault(sum_ln2piVariances) || numtype(sum_ln2piVariances)) //can provied sum_ln2pivariances externally, for speed
			MatrixOp/O/FREE ln2piVariances = ln(2*pi*varianceWave)
			Wavestats/Q/M=1 ln2piVariances
			sum_ln2piVariances = v_sum
		endif
		Npar -= 1 //variance is not a fit parameter
		SIC = (Npar)*ln(n) + sum_ln2piVariances + sum_scaledresidualssq
		
		Variable/G fld:noiseSigma = NaN
		Variable/G fld:noiseSigmaEstimated = NaN
		Variable/G fld:redchisq = 1/(n-Npar) * sum_scaledresidualssq
	endif
	
	return SIC
End



//Fitting with constraints is not available for the built-in line fit
Threadsafe Function LineFitOffset(Wave w, Variable x) : FitFunc
	//y = k*(x-x0)+y0
	return w[0]*(x-w[1])+w[2]
End


Threadsafe Static Function SortPauses(Wave w)
	if(dimsize(w,0)==0)
		return NaN
	endif
	
	Duplicate/O/FREE w, cpy
	Make/O/N=(dimsize(w,0))/FREE idx = p
	Duplicate/O/R=[][0] w, starts
	sort starts, starts, idx
	w[][] = cpy[idx[p]][q]
End


//Threadsafe Static Function getSlope(Wave w, Variable fromPnt, Variable toPnt [, Wave xwave])
//	Duplicate/O/FREE w, wX, wXY
//	
//	if(!paramisdefault(xwave) && waveExists(xwave))
//		wX = !numtype(w) ? xwave[p] : NaN
//	else
//		wX = !numtype(w) ? x : NaN
//	endif
//	
//	if(toPnt-fromPnt < 1)
//		return NaN
//	endif
//	
//	wXY = w*wX
//	Wavestats/Q/M=0/R=[fromPnt,toPnt] w
//	Variable SY = v_sum
//	
//	Wavestats/Q/M=0/R=[fromPnt,toPnt] wXY
//	Variable SXY = v_sum
//	
//	Wavestats/Q/M=2/R=[fromPnt,toPnt] wX
//	Variable SX = v_sum
//	Variable SX2 = v_rms^2*v_npnts
//	
//	return (SXY-SX*SY/v_npnts)/(SX2-SX^2/v_npnts)
//End


Threadsafe Static Function getSlope(Wave w, Variable fromPnt, Variable toPnt [, Wave xwave, Wave varwave])
    Duplicate/O/FREE w, wX, wXY
    Redimension/D wX, wXY
    Variable useWeights = !ParamIsDefault(varwave) && WaveExists(varwave)

    if(!ParamIsDefault(xwave) && WaveExists(xwave))
        wX = !numtype(w) ? xwave[p] : NaN
    else
        wX = !numtype(w) ? x : NaN
    endif

    if(toPnt - fromPnt < 1)
        return NaN
    endif

    if(useWeights)
        Duplicate/O/FREE w, wW, wWX, wWY, wWXY, wWX2
        Redimension/D wW, wWX, wWY, wWXY, wWX2
        wW   = (varwave[p] > 0 && !numtype(varwave[p])) ? 1/varwave[p] : NaN
        wWX  = wW * wX
        wWY  = wW * w
        wWXY = wW * wX * w
        wWX2 = wW * wX^2

        Wavestats/Q/M=0/R=[fromPnt,toPnt] wW;   Variable SW   = V_sum
        Wavestats/Q/M=0/R=[fromPnt,toPnt] wWX;  Variable SWX  = V_sum
        Wavestats/Q/M=0/R=[fromPnt,toPnt] wWY;  Variable SWY  = V_sum
        Wavestats/Q/M=0/R=[fromPnt,toPnt] wWXY; Variable SWXY = V_sum
        Wavestats/Q/M=0/R=[fromPnt,toPnt] wWX2; Variable SWX2 = V_sum

        return (SW*SWXY - SWX*SWY) / (SW*SWX2 - SWX^2)
    else
        wXY = w * wX
        Wavestats/Q/M=0/R=[fromPnt,toPnt] w;   Variable SY = V_sum
        Wavestats/Q/M=0/R=[fromPnt,toPnt] wXY; Variable SXY = V_sum
        Wavestats/Q/M=2/R=[fromPnt,toPnt] wX
        Variable SX = V_sum, SX2 = V_rms^2 * V_npnts

        return (SXY - SX*SY/V_npnts) / (SX2 - SX^2/V_npnts)
    endif
End


//Look-up table for end-slopes
Threadsafe Static Function MakeEndSlopeLUT(Wave w [, Wave xwave, Wave variancewave])
	DFREF fld = root:Packages:FindPauses
	
	Duplicate/O/FREE/D w, wX
	Variable useWeights = !ParamIsDefault(variancewave) && WaveExists(variancewave)
	Variable useXwave = !ParamIsDefault(xwave) && WaveExists(xwave)
	
	if(useXwave)
	    wX = !numtype(w) ? xwave[p] : NaN
	else
	    wX = !numtype(w) ? x : NaN
	endif

	Duplicate/O/FREE w, wW, wWX, wWY, wWXY, wWX2
	Redimension/D wW, wWX, wWY, wWXY, wWX2

	if(useWeights)
		wW  = !numtype(variancewave[p]) && !numtype(w) ? 1/variancewave[p] : 0
	else
		wW = !numtype(w) ? 1 : 0
	endif
	
	wWX  = !numtype(w) ? wW * wX : 0
	wWY  = !numtype(w) ? wW * w : 0
	wWXY = !numtype(w) ? wW * wX * w : 0
	wWX2 = !numtype(w) ? wW * wX^2 : 0
	
	Duplicate/O/FREE/D wW, SwW
	Duplicate/O/FREE/D wWX, SwWX
	Duplicate/O/FREE/D wWY, SwWY
	Duplicate/O/FREE/D wWXY, SwWXY
	Duplicate/O/FREE/D wWX2, SwWX2
	Integrate/P/METH=0 SwW, SwWX, SwWY, SwWXY, SwWX2
	
	Duplicate/O/D wW, fld:SlopeFromLeftEnd/Wave=SlopeFromLeftEnd
	SlopeFromLeftEnd = (SwW*SwWXY - SwWX*SwWY) / (SwW*SwWX2 - SwWX^2)
	
	Duplicate/O SwW, fld:SwW   //store to compute weighted means
	Duplicate/O SwWY, fld:SwWY //store to compute weighted means
	

//Duplicate/O/D wW, SlopeFromRightEnd
//slopefromrightend=nan
//Variable l = numpnts(w)-1
//SlopeFromRightEnd[1,] = ((SwW[l]-SwW[p-1])*(SwWXY[l]-SwWXY[p-1]) - (SwWX[l]-SwWX[p-1])*(SwWY[l]-SwWY[p-1])) / ((SwW[l]-SwW[p-1])*(SwWX2[l]-SwWX2[p-1]) - (SwWX[l]-SwWX[p-1])^2)

	//Because of numerical errors, I'll have to do it with reversal for computing from the other end. End-Integral gives slightly incorrect results just where I need accuracy.
	Duplicate/O/FREE/D wW, SwW
	Duplicate/O/FREE/D wWX, SwWX
	Duplicate/O/FREE/D wWY, SwWY
	Duplicate/O/FREE/D wWXY, SwWXY
	Duplicate/O/FREE/D wWX2, SwWX2
	Wavetransform	/O flip SwW
	Wavetransform	/O flip SwWX
	Wavetransform	/O flip SwWY
	Wavetransform	/O flip SwWXY
	Wavetransform	/O flip SwWX2
	Integrate/P/METH=0 SwW, SwWX, SwWY, SwWXY, SwWX2

	Duplicate/O/D wW, fld:SlopeFromRightEnd/Wave=SlopeFromRightEnd
	SlopeFromRightEnd = (SwW*SwWXY - SwWX*SwWY) / (SwW*SwWX2 - SwWX^2)
	Wavetransform	/O flip SlopeFromRightEnd
End


Threadsafe Static Function GetWeightedMean(Variable fromPnt, Variable toPnt)
	DFREF fld = root:Packages:FindPauses
	Wave/SDFR=fld SwW, SwWY
	Variable SWY = SwWY[toPnt] - (fromPnt > 0 ? SwWY[fromPnt-1] : 0)
	Variable SW  = SwW[toPnt]  - (fromPnt > 0 ? SwW[fromPnt-1]  : 0)
	return SWY / SW
End





///////////////////////////////// MCMC Version ///////////////////////////////////
//Do multiple runs with different random seeds

Function/WAVE FindPausesMCMC_MT(Wave traj [,Variable maxCountsSinceLastBest, Variable noiseSigma, Wave xwave, Variable randomseed, Wave sigmaWave, Variable Nruns, Variable quiet])
	maxCountsSinceLastBest = paramisdefault(maxCountsSinceLastBest) ? max(2000, dimsize(traj,0)*15) : maxCountsSinceLastBest //stop criterium. Bigger for better and slower.
	Nruns = paramisdefault(Nruns) ? 1 : Nruns
	randomseed = paramisdefault(randomseed) ? 0.1 : randomseed //default:reproducible
	quiet = paramisdefault(quiet) ? 0 : quiet
	
	if(Nruns<1)
		abort "Illegal Nruns"
	endif
	
	setrandomseed randomseed
	
	NewDataFolder/O root:packages
	NewDataFolder/O root:packages:FindPauses
	DFREF fld = root:packages:FindPauses
	
	Variable Nthreads = ThreadProcessorCount
	Nruns = min(Nthreads, Nruns)
	
	Variable threadGroupID = ThreadGroupCreate(Nthreads)
	Variable dummy, i
	for(i=0;i<Nruns;i+=1)
		//check if free thread index is available
		Variable threadIdx = ThreadGroupWait(threadGroupID,-2)-1
		if(threadIdx<0)
			dummy = ThreadGroupWait(threadGroupID, 50) //give threads some time
			i-=1
			continue
		endif
		
		//start thread
		ThreadStart threadGroupID, threadIdx, FindPausesMCMC_singlerun(traj, maxCountsSinceLastBest, noiseSigma, xwave, enoise(.5)+.5, sigmaWave, 1, 1)
	endfor
	
	//wait for all to finish
	try
		do
			DoUpdate
		while(ThreadGroupWait(threadGroupID, 100)!=0)
	catch
		dummy = ThreadGroupRelease(threadGroupID)
		abort
	endtry
	
	Make/O/FREE/N=(Nruns) SICs
	Make/O/FREE/DF/N=(Nruns) DFs

	//retrieve data via free data folder queue
	for(i=0;i<Nruns;i+=1)
		DFREF dfr= ThreadGroupGetDFR(threadGroupID,100)
				
		NVAR/SDFR=dfr SIC
		SICs[i] = SIC
		DFs[i] = dfr
	endfor	

	//get result with minimum SIC and store it in global data folder
	Wavestats/Q SICs
	DFREF bestdf = DFs[v_minloc]
	MoveDataFolder/O=3 bestdf, root:packages
	
	dummy = ThreadGroupRelease(threadGroupID)
	
	
	//Analyze result and print
	if(!paramisdefault(sigmawave) && WaveExists(sigmawave))
		//abort "Not implemented yet. This is more difficult than I thought. Need to fix fitting (weighted avg) and end-slopes, as well as maxLik."
		Duplicate/O/FREE sigmawave, variancewave
		variancewave = !numtype(traj) ? sigmawave^2 : NaN

		Duplicate/O sigmawave, fld:sigmawave
		MatrixOp/O/FREE ln2piVariances = ln(2*pi*varianceWave)
		Wavestats/Q/M=1 ln2piVariances
		Variable sum_ln2piVariances = v_sum
	else
		sum_ln2piVariances = NaN
	endif

	Wave/SDFR=fld W_Pauses
	Variable Npauses = dimsize(W_Pauses,0)
	SIC = GetSIC(traj, W_Pauses, noiseSigma=noiseSigma, xwave=xwave, variancewave=variancewave, sum_ln2piVariances=sum_ln2piVariances)
	NVAR NoiseSigmaEstimated = fld:noiseSigmaEstimated, NoiseSigmaOutput=fld:NoiseSigma, redChisq=fld:redChisq
	if(Npauses==1 && W_Pauses[0][0] == 0 && W_Pauses[0][1] == dimsize(traj,0)-1)
		String notice = "[Only one pause => Is static]"
		Variable/G fld:V_isStatic = 1
		Note W_Pauses, "isStatic:1"
	else
		Variable/G fld:v_isStatic = 0
		Note W_Pauses, "isStatic:0"
		notice = ""
	endif
	
	if(!quiet)
		printf "%d pause(s) found. SIC=%f, noiseSigmaEst=%f, noiseSigma=%f, redChisq=%f %s\r", Npauses, SIC, NoiseSigmaEstimated, NoiseSigmaOUtput, redchisq, notice
	endif
		
	return fld:W_pauses	
End




Function/Wave FindPausesMCMC(Wave traj [,Variable maxCountsSinceLastBest, Variable noiseSigma, Wave xwave, Variable randomseed, Wave sigmaWave, Variable Nruns, Variable quiet])
	maxCountsSinceLastBest = paramisdefault(maxCountsSinceLastBest) ? max(2000, dimsize(traj,0)*5) : maxCountsSinceLastBest //stop criterium. Bigger for better and slower.
	Nruns = paramisdefault(Nruns) ? 1 : Nruns
	randomseed = paramisdefault(randomseed) ? 0.1 : randomseed //default:reproducible
	quiet = paramisdefault(quiet) ? 0 : quiet
	
	if(Nruns<1)
		abort "Illegal Nruns"
	endif
	
	setrandomseed randomseed
	
	NewDataFolder/O root:packages
	NewDataFolder/O root:packages:FindPauses
	DFREF fld = root:packages:FindPauses
	
		
	Make/O/FREE/N=(Nruns) SICs
	
	try
		NewDataFolder root:TMP
		Variable i
		for(i=0;i<Nruns;i+=1)
			FindPausesMCMC_singlerun(traj, maxCountsSinceLastBest, noiseSigma, xwave, enoise(.5)+.5, sigmaWave, 1, 0)
			NVAR/SDFR=fld SIC
			
			DuplicateDataFolder fld, $("root:TMP:RUN_"+num2istr(i))
			
			SICs[i] = SIC
		endfor
		
		Wavestats/Q SICs
		DuplicateDataFolder/O=3 $("root:TMP:RUN_"+num2istr(v_minloc)), root:packages:FindPauses

		//clean up
		for(i=0;i<Nruns;i+=1)
			KillDataFolder/Z $("root:TMP:RUN_"+num2istr(i))
		endfor
		KillDataFolder TMP
	catch
		for(i=0;i<Nruns;i+=1)
			KillDataFolder/Z $("root:TMP:RUN_"+num2istr(i))
		endfor
		KillDataFolder root:TMP
	endtry

	if(!paramisdefault(sigmawave) && WaveExists(sigmawave))
		//abort "Not implemented yet. This is more difficult than I thought. Need to fix fitting (weighted avg) and end-slopes, as well as maxLik."
		Duplicate/O/FREE sigmawave, variancewave
		variancewave = !numtype(traj) ? sigmawave^2 : NaN

		Duplicate/O sigmawave, fld:sigmawave
		MatrixOp/O/FREE ln2piVariances = ln(2*pi*varianceWave)
		Wavestats/Q/M=1 ln2piVariances
		Variable sum_ln2piVariances = v_sum
	else
		sum_ln2piVariances = NaN
	endif

	Wave/SDFR=fld W_Pauses
	Variable Npauses = dimsize(W_Pauses,0)
	SIC = GetSIC(traj, W_Pauses, noiseSigma=noiseSigma, xwave=xwave, variancewave=variancewave, sum_ln2piVariances=sum_ln2piVariances)
	NVAR NoiseSigmaEstimated = fld:noiseSigmaEstimated, NoiseSigmaOutput=fld:NoiseSigma, redChisq=fld:redChisq
	if(Npauses==1 && W_Pauses[0][0] == 0 && W_Pauses[0][1] == dimsize(traj,0)-1)
		String notice = "[Only one pause => Is static]"
		Variable/G fld:V_isStatic = 1
		Note W_Pauses, "isStatic:1"
	else
		Variable/G fld:v_isStatic = 0
		Note W_Pauses, "isStatic:0"
		notice = ""
	endif
	
	if(!quiet)
		printf "%d pause(s) found. SIC=%f, noiseSigmaEst=%f, noiseSigma=%f, redChisq=%f %s\r", Npauses, SIC, NoiseSigmaEstimated, NoiseSigmaOUtput, redchisq, notice
	endif

	//re-seed RNG
	SetRandomSeed (ticks + trunc(1e6*abs(enoise(1))))/1e9
		
	return fld:W_pauses
End



//NOTES:
// - Zero-length pauses are fishy (change points). Unclear how to deal with them. Currently disallowing!!!. They will always sit exactly on data points (see GetSIC()). So this isn't a true change-point analysis.
// - xwave: optional x-wave for traj. If not given, wave scaling is used.
// - Fitting starts with no pause, then the first move is hard-coded to birth a pause that encompasses the entire traj. This way, the "static" and "moving" conditions are always tested.
Threadsafe Static Function FindPausesMCMC_singlerun(Wave traj, Variable maxCountsSinceLastBest, Variable noiseSigma, Wave xwave, Variable randomSeed, Wave sigmaWave, Variable quiet, Variable isThread)

	NewDataFolder/O root:packages
	NewDataFolder/O root:packages:FindPauses
	DFREF fld = root:packages:FindPauses
	
	Duplicate/O traj, fld:traj/Wave=w_traj
	
	if(WaveExists(xwave))
		Duplicate/O xwave, fld:xwave
	endif
	if(WaveExists(sigmawave))
		//abort "Not implemented yet. This is more difficult than I thought. Need to fix fitting (weighted avg) and end-slopes, as well as maxLik."
		Duplicate/O/FREE sigmawave, variancewave
		variancewave = !numtype(traj) ? sigmawave^2 : NaN

		Duplicate/O sigmawave, fld:sigmawave
		MatrixOp/O/FREE ln2piVariances = ln(2*pi*varianceWave)
		Wavestats/Q/M=1 ln2piVariances
		Variable sum_ln2piVariances = v_sum
	else
		sum_ln2piVariances = NaN
	endif
	MakeEndSlopeLUT(traj, xwave=xwave, variancewave=variancewave)
	
	SetRandomSeed randomseed
	
	Variable NMOVES = 6 //birth, death, changeLength, split, merge, fillToEnd, (changeLengthMinimal)
	
	Make/O/N=(NMOVES, 3) fld:P_MOVE/Wave=P_MOVE = 0 //second index: 0: no pauses, 1: one pause, 2: two pauses or more
	P_MOVE[][0] = {1,0,0,0,0,0}                     //zero pauses: only birth allowed
	P_MOVE[][1] = {1,1,1,.1,0,.05}                  //one pause: no merge
	P_MOVE[][2] = {.8, .2, 1, .1, .1, .02}          //two pauses or more: all allowed
	MatrixOp/O P_MOVE = P_MOVE / RowRepeat(sumcols(P_MOVE), NMOVES) //Normalize P_MOVE

	Duplicate/O P_MOVE, fld:CUM_P_MOVE/Wave = CUM_P_MOVE
	MatrixOp/O CUM_P_MOVE = integrate(P_MOVE, 1)

	Make/O/N=(0,2) fld:W_CurrentPauses/Wave=W_CurrentPauses
	SetDimLabel 1, 0, StartP, W_CurrentPauses
	SetDimLabel 1, 1, LengthP, W_CurrentPauses
	
	Make/O/N=0 fld:timeline/Wave=timeline
	Make/O/N=(NMOVES,6) fld:Movestats/Wave=Movestats = 0 //reject, accept
	SetDimLabel 0, 0, Birth, Movestats
	SetDimLabel 0, 1, Death, Movestats
	SetDimLabel 0, 2, ChangeLength, Movestats
	SetDimLabel 0, 3, Split, Movestats
	SetDimLabel 0, 4, Merge, Movestats
	SetDimlabel 0, 5, FillToEnd, Movestats
	//SetDimLabel 0, 6, ChangeLengthMinimal, Movestats
	SetDimLabel 1, 0, Accept, Movestats
	SetDimLabel 1, 1, Reject, Movestats
	SetDimLabel 1, 2, InfFW, Movestats
	SetDimLabel 1, 3, InfREV, Movestats
	SetDimLabel 1, 4, tComplete, Movestats //time for move completion (result not -inf)
	SetDimLabel 1, 5, tIncomplete, Movestats //time for move (result -inf)
	
	//---Initialize
	Duplicate/O W_CurrentPauses, fld:W_Pauses
	Variable lnp_w_given_Y = -GetSIC(W_traj, W_CurrentPauses, noiseSigma=noiseSigma, xwave=xwave, variancewave=variancewave, sum_ln2piVariances=sum_ln2piVariances)	
	InsertPoints/V=(lnp_w_given_Y) numpnts(timeline), 1, timeline
	Duplicate/O fld:fit, fld:bestFit/Wave=bestfit, fld:bestFit_colorwave/Wave=bestfit_colorwave, fld:bestfit_residuals/Wave=bestfit_residuals
	bestfit_colorwave = 0
	bestfit_residuals = W_traj-bestfit

	Variable ii
	for(ii=0;ii<dimsize(W_currentPauses,0);ii+=1)
		bestfit_colorwave[W_currentPauses[ii][0],max(0,W_currentPauses[ii][0]+W_currentPauses[ii][1]-1)] = 1
	endfor

	Variable lnp_wp_given_Y = NaN
	Variable lnProposal_w_given_wp = NaN
	Variable lnProposal_wp_given_w = NaN
	Variable best = lnp_w_given_Y
	Variable counterLastBest = -1
	Variable lastUpdate = ticks
	
	//---Loop
	Variable counter
	for(counter=0;;counter+=1)
		Variable move = SelectMove(CUM_P_MOVE, W_CurrentPauses)

		if(counter==0)
			move = -1 //special move first: full pause
		endif

		switch(move)
			case -1: //----- SPECIAL MOVE (only at counter=0) BIRTH OF FULL PAUSE ------
				//Generate proposal
				DFREF fld = root:Packages:FindPauses
				Duplicate/O W_CurrentPauses, fld:W_ProposedPauses/Wave=W_ProposedPauses
				if(dimsize(W_ProposedPauses,0)==0)
					make/O/N=(1,2) fld:W_ProposedPauses/Wave=W_ProposedPauses
					SetDimLabel 1, 0, StartP, W_ProposedPauses
					SetDimLabel 1, 1, LengthP, W_ProposedPauses
					W_ProposedPauses[0][0] = 0
					W_ProposedPauses[0][1] = dimsize(W_traj,0)-1	
				else
					print "Cannot do special move"
					return NaN
				endif
				
				lnProposal_wp_given_w = 0
				lnProposal_w_given_wp = 0
				break	
		
			case 0: //----- BIRTH ------
				lnProposal_wp_given_w = Move_Birth(w_traj, W_CurrentPauses)
				Wave/SDFR=fld W_ProposedPauses
				
				if(lnProposal_wp_given_w == -inf)
					Movestats[move][%InfFW] += 1
					break
				endif
	
				lnProposal_w_given_wp = lnP_death(w_traj, W_ProposedPauses)
	
				if(lnProposal_w_given_wp == -inf)
					Movestats[move][%InfRev] += 1
					break
				endif
	
				lnProposal_wp_given_w += ln(Get_P_Move(W_CurrentPauses, 0))
				lnProposal_w_given_wp += ln(Get_P_Move(W_ProposedPauses, 1))
				break


			case 1: //----- DEATH ------
				lnProposal_wp_given_w = Move_Death(w_traj, W_CurrentPauses)
				Wave/SDFR=fld W_ProposedPauses
				NVAR/SDFR=fld DEATHINDEX
				
				if(lnProposal_wp_given_w == -inf)
					Movestats[move][%InfFW] += 1
					break
				endif
	
				lnProposal_w_given_wp = lnP_birth(w_traj, W_CurrentPauses, W_ProposedPauses, DEATHINDEX)
	
				if(lnProposal_w_given_wp == -inf)
					Movestats[move][%InfRev] += 1
					break
				endif
	
				lnProposal_wp_given_w += ln(Get_P_Move(W_CurrentPauses, 1))
				lnProposal_w_given_wp += ln(Get_P_Move(W_ProposedPauses, 0))
				break

			case 2: //----- CHANGELENGTH ------
				lnProposal_wp_given_w = Move_ChangeLength(W_traj, W_CurrentPauses)
				Wave/SDFR=fld W_ProposedPauses
				NVAR/SDFR=fld CHANGELENGTH_INDEX, CHANGELENGTH_SIDE
				
				if(lnProposal_wp_given_w == -inf)
					Movestats[move][%InfFW] += 1
					break
				endif
	
				lnProposal_w_given_wp = lnP_changelength(W_traj, W_ProposedPauses, CHANGELENGTH_INDEX, CHANGELENGTH_SIDE)
	
				if(lnProposal_w_given_wp == -inf)
					Movestats[move][%InfRev] += 1
					break
				endif
	
				lnProposal_wp_given_w += ln(Get_P_Move(W_CurrentPauses, 2))
				lnProposal_w_given_wp += ln(Get_P_Move(W_ProposedPauses, 2))
				break

			case 3: //----- SPLIT ------
				lnProposal_wp_given_w = Move_Split(W_traj, W_CurrentPauses)
				Wave/SDFR=fld W_ProposedPauses
				
				if(lnProposal_wp_given_w == -inf)
					Movestats[move][%InfFW] += 1
					break
				endif
	
				lnProposal_w_given_wp = lnP_merge(W_traj, W_ProposedPauses)
	
				if(lnProposal_w_given_wp == -inf)
					Movestats[move][%InfRev] += 1
					break
				endif
	
				lnProposal_wp_given_w += ln(Get_P_Move(W_CurrentPauses, 3))
				lnProposal_w_given_wp += ln(Get_P_Move(W_ProposedPauses, 4))
				break

			case 4: //----- MERGE ------
				lnProposal_wp_given_w = Move_Merge(W_traj, W_CurrentPauses)
				Wave/SDFR=fld W_ProposedPauses
				NVAR/SDFR=fld MERGE_FIRSTLENGTH, MERGE_MERGEDLENGTH
				
				if(lnProposal_wp_given_w == -inf)
					Movestats[move][%InfFW] += 1
					break
				endif
	
				lnProposal_w_given_wp = lnP_split(W_traj, W_ProposedPauses, MERGE_MERGEDLENGTH, MERGE_FIRSTLENGTH)
	
				if(lnProposal_w_given_wp == -inf)
					Movestats[move][%InfRev] += 1
					break
				endif
	
				lnProposal_wp_given_w += ln(Get_P_Move(W_CurrentPauses, 4))
				lnProposal_w_given_wp += ln(Get_P_Move(W_ProposedPauses, 3))
				break

			case 5: //----- FILLTOEND ------
				lnProposal_wp_given_w = Move_FillToEnd(W_traj, W_CurrentPauses)
				Wave/SDFR=fld W_ProposedPauses
				NVAR/SDFR=fld FILLTOEND_INDEX, FILLTOEND_SIDE
				
				if(lnProposal_wp_given_w == -inf)
					Movestats[move][%InfFW] += 1
					break
				endif
	
				lnProposal_w_given_wp = lnP_changelength(W_traj, W_ProposedPauses, FILLTOEND_INDEX, FILLTOEND_SIDE)
	
				if(lnProposal_w_given_wp == -inf)
					Movestats[move][%InfRev] += 1
					break
				endif
	
				lnProposal_wp_given_w += ln(Get_P_Move(W_CurrentPauses, 5))
				lnProposal_w_given_wp += ln(Get_P_Move(W_ProposedPauses, 2))
				break

//			case 6: //----- CHANGELENGTHMINIMAL ------
//				lnProposal_wp_given_w = Move_ChangeLengthMinimal(W_traj, W_CurrentPauses)
//				Wave/SDFR=fld W_ProposedPauses
//				NVAR/SDFR=fld CHANGELENGTH_INDEX, CHANGELENGTH_SIDE
//				
//				if(lnProposal_wp_given_w == -inf)
//					Movestats[move][%InfFW] += 1
//					break
//				endif
//	
//				lnProposal_w_given_wp = lnP_changelengthMinimal(W_traj, W_ProposedPauses, CHANGELENGTH_INDEX, CHANGELENGTH_SIDE)
//	
//				if(lnProposal_w_given_wp == -inf)
//					Movestats[move][%InfRev] += 1
//					break
//				endif
//	
//				lnProposal_wp_given_w += ln(Get_P_Move(W_CurrentPauses, 6))
//				lnProposal_w_given_wp += ln(Get_P_Move(W_ProposedPauses, 6))
//				break

		endswitch
		
		//HACKY: move=-1 is special and will cause problems below. Set it to a "normal" birth move
		move = move==-1 ? 0 : move

		if(lnProposal_wp_given_w == -inf || lnProposal_w_given_wp == -inf)
			counter -= 1
			continue
		endif
		
		//Reject zero-length pauses
		if(dimsize(W_proposedPauses,0)>0)
			Wavestats/Q/RMD=[][1] W_proposedPauses
			if(v_min==0)
				counter -= 1
				continue
			endif
		endif

		lnp_wp_given_Y = -GetSIC(W_traj, W_ProposedPauses, noiseSigma=noiseSigma, xwave=xwave, variancewave=variancewave, sum_ln2piVariances=sum_ln2piVariances) //calculate posterior			

		Variable decision = MetropolisHastingsA(lnp_wp_given_Y, lnProposal_w_given_wp, lnp_w_given_Y, lnProposal_wp_given_w)

		if(decision==1)
			Duplicate/O W_ProposedPauses, fld:W_CurrentPauses
			lnp_w_given_Y = lnp_wp_given_Y //overwrite posterior
			Movestats[move][%Accept] += 1
		else
			Movestats[move][%Reject] += 1
		endif
		
		InsertPoints/V=(lnp_w_given_Y) numpnts(timeline), 1, timeline

		if(lnp_w_given_Y > best)
			best = lnp_w_given_Y
			counterlastbest = counter

			Wave/SDFR=fld fit
			Duplicate/O fit, fld:bestFit/Wave=bestfit
			Duplicate/O W_CurrentPauses, fld:W_Pauses
			
			//--Generate color wave (to mark pauses in bestfit) and residuals
			Duplicate/O fit, fld:bestfit_colorwave/Wave=bestfit_colorwave, fld:bestfit_residuals/Wave=bestfit_residuals
			bestfit_colorwave = 0
			bestfit_residuals = W_traj-bestfit
			for(ii=0;ii<dimsize(W_currentPauses,0);ii+=1)
				bestfit_colorwave[W_currentPauses[ii][0],max(0,W_currentPauses[ii][0]+W_currentPauses[ii][1]-1)] = 1
			endfor
			
		endif

		
		if(numtype(CheckPauses(W_traj, W_ProposedPauses, counter)) || numtype(CheckPauses(W_traj, w_currentPauses, counter)))
			return NaN
		endif
		
		Variable currentTicks = ticks
		if(0&&(currentTicks-lastUpdate)/60 > 1) //update every 1 secs?
			//doupdate //not available in Threadsafe functions
			lastUpdate = currentTicks
		endif
		
		if(counter-counterlastbest > maxCountsSinceLastBest)
			break
		endif
	endfor

	Wave/SDFR=fld W_Pauses
	Variable Npauses = dimsize(W_Pauses,0)
	Variable SIC = GetSIC(W_traj, W_Pauses, noiseSigma=noiseSigma, xwave=xwave, variancewave=variancewave, sum_ln2piVariances=sum_ln2piVariances)
	NVAR NoiseSigmaEstimated = fld:noiseSigmaEstimated, NoiseSigmaOutput=fld:NoiseSigma, redChisq=fld:redChisq
	if(Npauses==1 && W_Pauses[0][0] == 0 && W_Pauses[0][1] == dimsize(W_traj,0)-1)
		String notice = "[Only one pause => Is static]"
		Variable/G fld:V_isStatic = 1
		Note W_Pauses, "isStatic:1"
	else
		Variable/G fld:v_isStatic = 0
		Note W_Pauses, "isStatic:0"
		notice = ""
	endif
	
	if(!quiet)
		printf "%d pause(s) found. SIC=%f, noiseSigmaEst=%f, noiseSigma=%f, redChisq=%f %s\r", Npauses, SIC, NoiseSigmaEstimated, NoiseSigmaOUtput, redchisq, notice
	endif
	
	Variable/G fld:SIC = SIC
	
	//re-seed RNG
	SetRandomSeed (ticks + trunc(1e6*abs(enoise(1))))/1e9

	if(isThread)
		Waveclear W_pauses, bestfit_residuals, bestfit_colorwave, bestfit, fit, w_proposedPauses, w_currentPauses, CUM_P_MOVE, P_MOVE, w_traj, timeline, movestats //clear all references for ThreadGroupPutDF
		ThreadGroupPutDF 0, fld
	endif
	
	return SIC
End



Threadsafe Static Function Move_Birth(Wave W_traj, Wave W_CurrentPauses)
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses==0)
		Make/O/FREE NinPauses = {0}
	else
		MatrixOp/O/FREE NInPauses = sum(col(W_CurrentPauses,1)) + Npauses
	endif
	
	//Pick random pause start
	if(N-NinPauses[0]-1<=0) //no option any more?
		return -inf
	endif
	Variable Proposed_PauseStartIndex = IntNoise(0, N-NinPauses[0]-1)
	
	//Translate this index to traj-point-number
	Variable i, offset = 0
	for(i=0;i<Npauses;i+=1)
		if(Proposed_PauseStartIndex+offset < W_CurrentPauses[i][0])
			break
		endif
		offset += W_CurrentPauses[i][1]+1
	endfor //after this, i contains pause index before which new pause is to be inserted
	Variable Proposed_PauseStartPoint = Proposed_PauseStartIndex+offset

	//Pick random pause length
	Variable StartOfNextPause = i<Npauses ? W_CurrentPauses[i][0] : N
	Variable Proposed_PauseLength = IntNoise(0, StartOfNextPause-Proposed_PauseStartPoint-1)
	
	if(Proposed_PauseLength < 0)
		return -inf
	endif
	
	//Generate proposal
	DFREF fld = root:Packages:FindPauses
	Duplicate/O W_CurrentPauses, fld:W_ProposedPauses/Wave=W_ProposedPauses
	if(dimsize(W_ProposedPauses,0)==0)
		make/O/N=(0,2) fld:W_ProposedPauses/Wave=W_ProposedPauses
		SetDimLabel 1, 0, StartP, W_ProposedPauses
		SetDimLabel 1, 1, LengthP, W_ProposedPauses
	endif

	InsertPoints/M=0 i, 1, W_ProposedPauses
	W_ProposedPauses[i][0] = Proposed_PauseStartPoint
	W_ProposedPauses[i][1] = Proposed_PauseLength	
	
	return -ln(N-NinPauses[0]) - ln(StartOfNextPause-Proposed_PauseStartPoint-1) 
end



Threadsafe Static Function LnP_Birth(Wave W_traj, Wave W_PausesBeforeDeath, Wave W_PausesAfterDeath, Variable KilledPauseNumber)
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_PausesAfterDeath,0)
	if(Npauses==0)
		Make/O/FREE NinPauses = {0}
	else
		MatrixOp/O/FREE NInPauses = sum(col(W_PausesAfterDeath,1)) + Npauses
	endif

	Variable KilledPause_Start = W_PausesBeforeDeath[KilledPauseNumber][0]
	Variable KilledPause_Length = W_PausesBeforeDeath[KilledPauseNumber][1]

	Variable NextPause_Start = KilledPauseNumber<dimsize(W_PausesAfterDeath,0) ? W_PausesAfterDeath[KilledPauseNumber][0] : N
	
	return -ln(N-NinPauses[0]) -ln(NextPause_Start-KilledPause_Start-1)
End



Threadsafe Static Function Move_Death(Wave W_traj, Wave W_CurrentPauses)
	Variable rndPause = IntNoise(0,dimsize(W_CurrentPauses,0)-1)

	if(dimsize(W_CurrentPauses,0)==0)
		return -inf
	endif

	DFREF fld = root:Packages:FindPauses
	Duplicate/O W_CurrentPauses, fld:W_ProposedPauses/Wave=W_ProposedPauses
	DeletePoints/M=0 rndPause, 1, W_ProposedPauses

	Variable/G fld:DEATHINDEX = rndPause

	return -ln(dimsize(W_CurrentPauses,0))
End


Threadsafe Static Function LnP_Death(Wave traj, Wave W_CurrentPauses)
	return -ln(dimsize(W_CurrentPauses,0))
End




Threadsafe Static Function Move_ChangeLength(Wave W_traj, Wave W_CurrentPauses)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses==0)
		return -inf
	endif
	
	//--Pick random pause
	Variable index = IntNoise(0, Npauses-1)
	Variable/G fld:CHANGELENGTH_INDEX = index
	lnP -= ln(Npauses)
	
	//--Pick side (+1: right side, -1: left side)
	Variable side = enoise(1)>0 ? 1 : -1
	Variable/G fld:CHANGELENGTH_SIDE = side
	lnP -= ln(2)
	
	//--Generate Proposal
	Duplicate/O W_CurrentPauses, fld:W_ProposedPauses/Wave=W_ProposedPauses
	if(side == +1) //Change right end?
		Variable startOfThisPause = W_CurrentPauses[index][0]
		Variable startOfNextPause = index<Npauses-1 ? W_CurrentPauses[index+1][0] : N
		
		Variable newLength = IntNoise(0, startOfNextPause-startOfThisPause-1)
		lnP -= ln(startOfNextPause-startOfThisPause-1)
		
		W_ProposedPauses[index][1] = newLength
	else //Change of left end?
		Variable endOfThisPause = W_CurrentPauses[index][0] + W_CurrentPauses[index][1]
		Variable endOfPreviousPause = index>0 ? W_CurrentPauses[index-1][0]+W_CurrentPauses[index-1][1] : -1
		
		Variable newStart = IntNoise(endOfPreviousPause+1, endOfThisPause)
		lnP -= ln(endOfThisPause-(endOfPreviousPause+1))
		newLength = endOfThisPause-newStart
		
		W_ProposedPauses[index][0] = newStart
		W_ProposedPauses[index][1] = newLength
	endif
	
	return lnP
end



Threadsafe Static Function LnP_ChangeLength(Wave W_traj, Wave W_CurrentPauses, Variable index, Variable side)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses==0)
		return -inf
	endif
	
	//--Pick random pause
	lnP -= ln(Npauses)
	
	//--Pick direction
	lnP -= ln(2)
	
	//--Generate Proposal
	if(side == +1) //Change right end?
		Variable startOfThisPause = W_CurrentPauses[index][0]
		Variable startOfNextPause = index<Npauses-1 ? W_CurrentPauses[index+1][0] : N
		
		lnP -= ln(startOfNextPause-startOfThisPause-1)
	else //Change of left end?
		Variable endOfThisPause = W_CurrentPauses[index][0] + W_CurrentPauses[index][1]
		Variable endOfPreviousPause = index>0 ? W_CurrentPauses[index-1][0]+W_CurrentPauses[index-1][1] : -1
		
		Variable newStart = IntNoise(endOfPreviousPause+1, endOfThisPause)
		lnP -= ln(endOfThisPause-(endOfPreviousPause+1))
	endif
	
	return lnP
end



//Static Constant CHANGELENGTH_CHARACTERISTICLENGTH = 3 //points

Threadsafe Static Function Move_ChangeLengthMinimal(Wave W_traj, Wave W_CurrentPauses)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses==0)
		return -inf
	endif
	
	//--Pick random pause
	Variable index = IntNoise(0, Npauses-1)
	Variable/G fld:CHANGELENGTH_INDEX = index
	lnP -= ln(Npauses)
	
	//--Pick side (+1: right side, -1: left side)
	Variable side = enoise(1)>0 ? 1 : -1
	Variable/G fld:CHANGELENGTH_SIDE = side
	lnP -= ln(2)

	//--Pick extension/shrinkage length
	Variable extendOrShrinkLength = enoise(1)>0 ? 1 : -1	

	Variable startOfThisPause = W_CurrentPauses[index][0]
	Variable endOfThisPause = W_CurrentPauses[index][0] + W_CurrentPauses[index][1]
	
	//--Generate Proposal
	Duplicate/O W_CurrentPauses, fld:W_ProposedPauses/Wave=W_ProposedPauses
	if(side == +1) //Change right end?
		Variable startOfNextPause = index<Npauses-1 ? W_CurrentPauses[index+1][0] : N
		
		Variable newLength = endOfThisPause-startOfThisPause + extendOrShrinkLength
		
		if(newLength<0 || startOfThisPause+newLength >= startOfNextPause)
			return -inf
		endif
		
		W_ProposedPauses[index][1] = newLength
	else //Change of left end?
		Variable endOfPreviousPause = index>0 ? W_CurrentPauses[index-1][0]+W_CurrentPauses[index-1][1] : -1
		
		Variable newStart = startOfThisPause + extendOrShrinkLength
		newLength = endOfThisPause-newStart

		if(newStart<0 || newStart <= endOfPreviousPause)
			return -inf
		endif
		
		W_ProposedPauses[index][0] = newStart
		W_ProposedPauses[index][1] = newLength
	endif
	
	return lnP
end

Threadsafe Static Function LnP_ChangeLengthMinimal(Wave W_traj, Wave W_CurrentPauses, Variable index, Variable side)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses==0)
		return -inf
	endif
	
	//--Pick random pause
	lnP -= ln(Npauses)
	
	//--Pick direction
	lnP -= ln(2)
	
	//--Pick extension/shrinkage length
	lnP -= ln(2)
	
	return lnP
end






//randomly select first or last pause, then make pause go from the start or to the end.
Threadsafe Static Function Move_FillToEnd(Wave W_traj, Wave W_CurrentPauses)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses==0)
		return -inf
	endif
	
	//--Pick random pause
	if(Npauses==1)
		Variable index = 0 //only one: pick it
		lnP -= 0
		Variable/G fld:FILLTOEND_NSELECTFROM = 1
	else
		index = enoise(1)>0 ? 0 : Npauses-1 //pick either first or last
		lnP -= ln(2)
		Variable/G fld:FILLTOEND_NSELECTFROM = 2
	endif	
	Variable/G fld:FILLTOEND_INDEX = index
	
	//--Pick side (+1: right side, -1: left side)
	if(Npauses==1)
		Variable side = enoise(1)>0 ? 1 : -1
		Variable/G fld:FILLTOEND_SIDE = side
		lnP -= ln(2)
	else
		if(index==0)
			side = -1
			Variable/G fld:FILLTOEND_SIDE = -1
		elseif(index==Npauses-1)
			side = +1
			Variable/G fld:FILLTOEND_SIDE = +1
		else
			print "ERROR"
			return NaN
		endif
	endif
	
	//--Generate Proposal
	Duplicate/O W_CurrentPauses, fld:W_ProposedPauses/Wave=W_ProposedPauses
	if(side == +1) //Change right end?
		Variable startOfThisPause = W_CurrentPauses[index][0]
		Variable newLength = N-startOfThisPause-1
		
		W_ProposedPauses[index][1] = newLength
	else //Change of left end?
		Variable endOfThisPause = W_CurrentPauses[index][0] + W_CurrentPauses[index][1]
		
		Variable newStart = 0
		newLength = endOfThisPause-newStart
		
		W_ProposedPauses[index][0] = newStart
		W_ProposedPauses[index][1] = newLength
	endif
	
	return lnP
end



Threadsafe Static Function Move_Merge(Wave W_traj, Wave W_CurrentPauses)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses<2)
		return -inf
	endif
	
	//--Pick random pause (first pause)
	Variable indexLeft = IntNoise(0, Npauses-2)
	Variable/G fld:MERGE_INDEXLEFT = indexLeft
	lnP -= ln(Npauses-1)
		
	//--Generate Proposal
	Duplicate/O W_CurrentPauses, fld:W_ProposedPauses/Wave=W_ProposedPauses
	Variable startOfLeftPause = W_CurrentPauses[indexLeft][0]
	Variable endOfRightPause = W_CurrentPauses[indexLeft+1][0] + W_CurrentPauses[indexLeft+1][1]
	Variable mergedLength = endOfRightPause - startOfLeftPause

	Variable/G fld:MERGE_FIRSTLENGTH = W_CurrentPauses[indexLeft][1], fld:MERGE_MERGEDLENGTH = mergedLength

	DeletePoints/M=0 indexLeft+1, 1, W_ProposedPauses
	W_ProposedPauses[indexLeft][0] = startOfLeftPause
	W_ProposedPauses[indexLeft][1] = mergedLength
	
	return lnP
end


Threadsafe Static Function LnP_Merge(Wave W_traj, Wave W_CurrentPauses)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses<2)
		return -inf
	endif
	
	lnP -= ln(Npauses-1)
	
	return lnP
end



Threadsafe Static Function Move_Split(Wave W_traj, Wave W_CurrentPauses)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses==0)
		return -inf
	endif
	
	//--Pick random pause
	Variable index = IntNoise(0, Npauses-1)
	lnP -= ln(Npauses)
	
	Variable CurrentLength = W_CurrentPauses[index][1]
	if(CurrentLength==0)
		return -inf
	endif

	//--Generate Proposal
	Duplicate/O W_CurrentPauses, fld:W_ProposedPauses/Wave=W_ProposedPauses
	
	//Pick new length first pause
	Variable firstStart = W_CurrentPauses[index][0]
	Variable firstLength = IntNoise(0, CurrentLength-1)
	lnP -= ln(CurrentLength)
	Variable/G fld:SPLIT_FIRSTLENGTH = firstLength
	
	//Pick new start of second pause
	Variable secondStart = IntNoise(W_CurrentPauses[index][0]+firstLength+1, W_CurrentPauses[index][0]+CurrentLength)
	lnP -= ln(CurrentLength-firstLength)
	Variable mergedLength = firstStart+CurrentLength-secondStart

	W_ProposedPauses[index][0] = firstStart //stays the same
	W_ProposedPauses[index][1] = firstLength
	
	InsertPoints/M=0 index+1, 1, W_ProposedPauses
	W_ProposedPauses[index+1][0] = secondStart
	W_ProposedPauses[index+1][1] = mergedLength
	
	return lnP
end


Threadsafe Static Function LnP_Split(Wave W_traj, Wave W_CurrentPauses, Variable LengthOfMergedPause, Variable LengthOfFirstPause)
	DFREF fld = root:Packages:FindPauses
	
	Variable lnP = 0
	
	Variable N = numpnts(W_traj)
	Variable Npauses = Dimsize(W_CurrentPauses,0)
	if(Npauses==0)
		return -inf
	endif
	
	//--Pick random pause
	lnP -= ln(Npauses)
	
	if(LengthOfMergedPause==0)
		return -inf
	endif

	//Pick new length first pause
	lnP -= ln(LengthOfMergedPause)
	
	//Pick new start of second pause
	lnP -= ln(LengthOfMergedPause-LengthOfFirstPause)
	
	return lnP
end



Threadsafe Static Function SelectMove(Wave CUM_P_MOVE, Wave W_Pauses)
	Variable Npauses = dimsize(W_Pauses,0)

	if(Npauses==0)
		Duplicate/O/FREE/R=[][0] CUM_P_MOVE, select
		Redimension/N=(-1) select
	elseif(Npauses==1)
		Duplicate/O/FREE/R=[][1] CUM_P_MOVE, select
		Redimension/N=(-1) select
	else
		Duplicate/O/FREE/R=[][2] CUM_P_MOVE, select
		Redimension/N=(-1) select
	endif

	return RandomIndexFromCDF(select)
End

Threadsafe Static Function Get_P_Move(Wave W_Pauses, Variable Move)
	DFREF fld = root:Packages:FindPauses
	Wave/SDFR=fld P_MOVE

	Variable Npauses = dimsize(W_Pauses,1)

	if(Npauses==0)
		return P_MOVE[Move][0]
	elseif(Npauses==1)
		return P_MOVE[Move][1]
	else
		return P_MOVE[Move][2]
	endif
End



Threadsafe Function CheckPauses(Wave W_traj, Wave W_Pauses, Variable counter)
	Variable Npauses = dimsize(W_Pauses,0), N = numpnts(W_traj)
	Variable i
	for(i=0;i<Npauses;i+=1)
		Variable startP = W_Pauses[i][0]
		Variable lengthP = W_Pauses[i][1]
		
		if(numtype(startP) || numtype(lengthP))
			print "CheckPauses failed. Numtype. " + num2istr(counter)
			return NaN
		endif
		if(startP<0 || lengthP<0 || startP >= N || startP+lengthP >= N)
			print "CheckPauses failed. Exceeded range. " + num2istr(counter)
			return NaN
		endif
		if(i<Npauses-1 && startP > w_Pauses[i+1][0])
			print "CheckPauses failed: Out of order. " + num2istr(counter)
			return NaN
		endif
		if(i<Npauses-1 && startP+lengthP >= w_Pauses[i+1][0])
			print "CheckPauses failed: Collision or overlap. " + num2istr(counter)
			return NaN
		endif
	endfor
	
	return 1
End

Threadsafe Static Function MetropolisHastingsA(Variable lnp_wp_given_Y, Variable lnProposal_w_given_wp, Variable lnp_w_given_Y, Variable lnProposal_wp_given_w)
	Variable lnNum   = lnp_wp_given_Y + lnProposal_w_given_wp
	Variable lnDenom = lnp_w_given_Y  + lnProposal_wp_given_w

	Variable A = min(1,exp(lnNum-lnDenom))
	Variable/G V_metropolisHastingsA = A
	Variable rnd = enoise(.5, 2)+.5
	if(rnd<A)
		return 1
	else
		return 0
	endif
End

// Generate uniformly-distributed integers on the interval [from,to] with from<to
Threadsafe Static Function IntNoise(from, to)
	Variable from, to
	Variable amp = to - from
	return floor(from + mod(abs(enoise(100*amp,2)),amp+1))
End

Threadsafe Static Function RandomIndexFromCDF(Wave CDF) //randomly pick from CDF (must be normalized!)
	Duplicate/O/FREE CDF, select

	if(abs(CDF[dimsize(CDF,0)-1] -1 ) > 1e-6)
		print "CDF not normalized"
		return NaN
		//abort "CDF not normalized"
	endif

	Variable rnd = enoise(.5, 2)+.5
	select = rnd<select
	Wavestats/Q/M=0 select
	return v_maxloc //ATTENTION: THIS ACTUALLY RETURNS THE SCALED POSITION AND NOT THE INDEX
End
