# IgorProCode

#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


////// Most important Igor Pro code written by June Groothuizen during bachelor project Physics and Astronomy 2020 /////


// finds an MDC at certain energy (input EDM)
Function getMDC(w,E)
	Wave w
	Variable E
	Variable Eindex = scaletoindex(w,E,1)
	
	Make/N=(dimsize(w,0))/O MDC
	SetScale/P x, dimoffset(w,0), dimdelta(w,0), MDC
	MDC = w[p][Eindex]
	
End


// Used to fit a Lorentian over an MDC
Function fitit(tofit)
	Wave tofit
	
	SetDataFolder root:global
	NVAR a, b, c, d, e, f, g, h, i	
	SetDataFolder root:fitting:result
	
	Make/O/N=(dimsize(tofit,0)) testshape
	SetScale/P x, dimoffset(tofit,0), dimdelta(tofit,0), testshape
	testshape = (1/(2 *pi)) * (a * (b/2) /((x-c)*(x-c) + (b/2)*(b/2))) +(1/(2 *pi)) * (d * (e/2) /((x-f)*(x-f) + (e/2)*(e/2) )) +g + x*h + x*x*i
	
	Make/O/N=9 w = {a, b, c, d, e, f, g, h, i}
	
	FuncFit/Q/H="000000111" fitfunction, w, tofit
End


// Fit for summation of 2 Lorentzians
Function fitfunction(w,x) : FitFunc
	Wave w
	Variable x
	
	Variable a = w[0]
	Variable b = w[1]
	Variable c = w[2]
	Variable d = w[3]
	Variable e = w[4]
	Variable f = w[5]
	Variable g = w[6]
	Variable h = w[7]
	Variable i = w[8]		
		
	Variable func = (1/(pi)) * (a * (b/2) /((x-c)*(x-c) + (b/2)*(b/2))) +(1/(pi)) * (d * (e/2) /((x-f)*(x-f) + (e/2)*(e/2) )) +g + x*h + x*x*i
	
	return func
End


// Function for fitting the MDC's for a whole EDM
function fitloop(EDM)

	Wave EDM
	Wave parameters
		
						// index		// energy
	Variable Ecut, E_final = 0, E_start = -0.003
	
	
	NewDatafolder/O/S root:fitting:result
	
	Make/O/N=(9,dimsize(EDM,1)) parameterwave
	Setscale/P y, dimoffset(EDM,1), dimdelta(EDM,1), parameterwave
	parameterwave = 0
	
	for (Ecut = scaletoindex(EDM,E_start,1); Ecut>(E_final-1); Ecut--)

		getMDC(EDM,indextoscale(EDM,Ecut,1))
		Wave MDC

		// determines the starting parameters using the MDC at 0
		if (Ecut == scaletoindex(EDM,E_start,1))
		
			Variable posL, peakL, HM_L, FWHM_L, width_L
			Wavestats/q/R=(-0.6,0) MDC
			posL = V_maxloc 
			peakL = V_max
			HM_L = V_max/2
			
			findvalue/V=(HM_L)/T=(0.09*HM_L) MDC
			width_L = indextoscale(MDC,V_value,0)
			FWHM_L = 2* abs(abs(posL) - abs(width_L))
			
			Variable posR, peakR, HM_R, FWHM_R, width_R
			Wavestats/q/R=(0,0.6) MDC
			posR = V_maxloc
			peakR = V_max
			HM_R = V_max/2

			findvalue/V=(HM_R)/T=(0.09*HM_R)/S=(scaletoindex(MDC,posR,0)) MDC
			width_R = indextoscale(MDC,V_value,0)
			FWHM_R = 2* abs(abs(posR) - abs(width_R))
			
			NewDataFolder/O/S root:global
			Variable/G a = peakL, b = FWHM_L, c = posL, d = peakR, e = FWHM_R, f = posR, g= 0, h = 0, i = 0	
			SetDataFolder root:Fitting:result
			
		endif
		
		// fits a Lorentizian over the MDC
		fitit(MDC)
		wave w
		
		// Updates the parameters for the next MDC
		SetDataFolder root:global
		Variable/G a = w[0], b = w[1], c = w[2], d = w[3], e = w[4], f = w[5], g= w[6], h = w[7], i = w[8]	
		SetDataFolder root:Fitting:result
		
		// fills the empty 2D wave with fitting values for each cut
		variable n
		
		for (n =0; n<9; n++)
			parameterwave[n][Ecut] = w[n]
		endfor
		
	endfor
	
End



// Plots values from the 2D parameterwave coming from the EDM fit
Function useparameterwave(vF,kf)
	Variable Vf, kf
	vf = -vf
	SetDataFolder root:fitting:result
	
	Wave parameterwave
	
	
	// Left part
	Make/O/N=(dimsize(parameterwave,1)) WL
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), WL
	WL = parameterwave[0][p]
	
	Make/O/N=(dimsize(parameterwave,1)) ImL
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), ImL
	//ImL = parameterwave[1][p] / (2 *parameterwave[0][p])  // with width W_L
	ImL = (parameterwave[1][p] * vF) / 2 // With some constant VF
	
	//ImL = parameterwave[1][p]
	
	Make/O/N=(dimsize(parameterwave,1)) ReL
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), ReL
	//ReL = - abs(kf - parameterwave[2][p])/ (parameterwave[0][p])  // With width
	//ReL = - abs(kf - parameterwave[2][p]) * 0.2 // With band assumptions
	ReL = parameterwave[2][p] //bare parameter
	
	// Right part
	Make/O/N=(dimsize(parameterwave,1)) WR
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), WR
	WR = parameterwave[3][p]
	
	Make/O/N=(dimsize(parameterwave,1)) ImR
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), ImR
	//ImR = parameterwave[4][p]/ (2 *parameterwave[3][p]) // width W_R
	ImR = (parameterwave[4][p] * vF) /2 // with some constant vF
	//ImR = parameterwave[4][p]
	
	Make/O/N=(dimsize(parameterwave,1)) ReR
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), ReR
	//ReR = - abs(kf - parameterwave[5][p])/ (parameterwave[3][p])  //With width
	ReR = parameterwave[5][p] // With band assumptions
	//ReR = parameterwave[5][p] //bare parameter
	
	// Background
	Make/O/N=(dimsize(parameterwave,1)) back_a
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), back_a
	back_a = parameterwave[6][p]
	Make/O/N=(dimsize(parameterwave,1)) back_b
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), back_b
	back_b = parameterwave[7][p]
	Make/O/N=(dimsize(parameterwave,1)) back_c
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), back_c
	back_c = parameterwave[8][p]

	//ModifyGraph rgb(ImR)=(13102,26214,0)
	//Display WL, ImL, ReL, WR, ImR, ReR, back_a, back_b, back_c
	//TextBox/C/N=Title/F=0/A=MC "Fitting result"
	//TextBox/C/N=Title/A=MT/X=0.00/Y=0.00/E
	//TextBox/C/N=text1/A=LT "IM:0.1, RE:-0.1, linear"
	//Label bottom "omega"
	//Label left "Intensity (a.u.)"
	//SetAxis left -0.4,0.6;DelayUpdate
	//SetAxis bottom -0.5,0.1
	//ModifyGraph rgb(ImR)=(19675,39321,1)
	//colourise(99,0)
	//legend
	
	// Can be used to check the values of a certain slice: plots the Lorentzian
	variable m = 208
	Make/O/N=(685) check
	Setscale/I x, -0.64, 0.76, check
	check = (1/ (pi)) * (parameterwave[0][m] * (parameterwave[1][m]/2) / ( (x-parameterwave[2][m])*(x-parameterwave[2][m]) + (parameterwave[1][m]/2)*(parameterwave[1][m]/2) )) + (1/ (pi)) * (parameterwave[3][m] * (parameterwave[4][m]/2) / ((x-parameterwave[5][m])*(x-parameterwave[5][m]) + (parameterwave[4][m]/2)*(parameterwave[4][m]/2) )) + parameterwave[6][m] + x * parameterwave[7][m] + x*x*parameterwave[8][m]

End


// Will make an EDM from a parameterwave as sanity check
Function checkparameterwave(parameterwave, EDM, REdisp, vf,kf)
	Wave parameterwave
	Wave EDM, REdisp
	Variable vf, kf
	
	NewDatafolder/O/S root:testfitting
	
	Duplicate/O/R=[,scaletoindex(EDM, 0,0)] EDM, EDM_split
	Duplicate/O EDM_split, EDM_test
	Duplicate/O EDM_split, EDM_values
	EDM_test = 0
	EDM_values = 0
	
	//Make/O/N=(100) check
	//Setscale/I x, -2, 2, check
	//check = (1/ (2 *pi)) * (parameterwave[0][m] * (parameterwave[1][m]/2) / ( (x-parameterwave[2][m])*(x-parameterwave[2][m]) + (parameterwave[1][m]/2)*(parameterwave[1][m]/2) )) + (1/ (2 *pi)) * (parameterwave[3][m] * (parameterwave[4][m]/2) / ((x-parameterwave[5][m])*(x-parameterwave[5][m]) + (parameterwave[4][m]/2)*(parameterwave[4][m]/2) )) + parameterwave[6][m] + x * parameterwave[7][m] + x*x*parameterwave[8][m]
	
	variable n
	//variable nfinal = scaletoindex(EDM,-1.2,1)
	variable nfinal = 0

	print scaletoindex(EDM_split,-0,1)-1
	for (n=(scaletoindex(EDM_split,-0,1)-1);n>nfinal;n=n-1)

		EDM_test[][n] = (1/ pi) * (parameterwave[0][n] * (parameterwave[1][n]/2) / ( (x-parameterwave[2][n])*(x-parameterwave[2][n]) + (parameterwave[1][n]/2)*(parameterwave[1][n]/2) )) + (1/ pi) * (parameterwave[3][n] * (parameterwave[4][n]/2) / ((x-parameterwave[5][n])*(x-parameterwave[5][n]) + (parameterwave[4][n]/2)*(parameterwave[4][n]/2) )) + parameterwave[6][n] + x * parameterwave[7][n] + x*x*parameterwave[8][n]
		
		EDM_values[][n] = (1/ pi) * (parameterwave[0][n] * vf * ((parameterwave[1][n]*vf)/2) / ( (indextoscale(EDM_split,n,1) - (vf *( x -kf)) - REdisp[dimsize(EDM_split,1)-1 - n])^2 + ((parameterwave[1][n]*vf)/2)^2 ) )
		
	endfor
	
	//Make/O/N=(100) check
	//Setscale/I x, -2, 2, check
	//check = 
	
	//Quantify: residual
	Duplicate/O EDM_split, residu
	residu = 0
	variable m, mfinal = 0
	for (m = scaletoindex(EDM,0,1);m>mfinal+1; m = m-1)
		//residu[][m] = EDM[p][m] - EDM_test[p][m]
		residu[][m] = EDM_split[p][m] - EDM_values[p][m]
		//F5_295K_GX = root:Fitting:F5_295K_GX[p][q] - (root:Fitting:F5_295K_GX_Hprofile0[p]/25)
	endfor
End


// Makes the dispersion for given bare band parameters
Function bare_dispersion(EDM,vf, kf)
	Variable vf, kf
	wave EDM
	
	SetDataFolder root:Fitting:KK
	wave rawfitpeakreversed, rawfitpeak
	
	variable last = indextoscale(EDM,0,1)/ vf + kf
	variable delta = -(kf - last) / dimsize(rawfitpeak,0)

	Make/O/N=(dimsize(rawfitpeak,0)) Dispersion
	SetScale/P x, kf, delta, Dispersion
	dispersion = vf * (x - kf)
End


// Real part from fit
Function rawRedata(parameterwave, EDM)
	Wave parameterwave
	Wave EDM
	
	SetDataFolder root:fitting:KK
	
	Make/O/N=(dimsize(parameterwave,1)) rawfitpeak_help
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), rawfitpeak_help
	
	
	Variable from0 = scaletoindex(EDM,0,1) -1

	rawfitpeak_help = parameterwave[2][p] //bare parameter
	Duplicate/O/R=[,from0] rawfitpeak_help, rawfitpeak
End


// Calculates real part of self-energy
Function ydifference(EDM)
	Wave EDM
	SetDataFolder root:Fitting:KK
	Wave dispersion, rawfitpeak, rawfitpeakreversed
	variable m = dimsize(dispersion,0), here = 300
	
	Variable ending =  m -scaletoindex(rawfitpeak, -0.6,0)
	
	Make/O/N=(ending) REdisp
	Setscale/P x, 0, dimdelta(EDM,1), REdisp
	variable n, fitE, fitk, dispind, dispE
	for (n= 0; n< ending; n++)
		fitE = indextoscale(rawfitpeak, m-n,0)
		fitk = rawfitpeak[m-n]
		dispind = scaletoindex(dispersion, fitK,0)
		dispE = dispersion[dispind]

		REdisp[n] = (abs(dispE) - abs(fitE))
	endfor
END


// Mirrors the IM part of the self-energy in E_F 
Function mirrorIM(w)
	Wave w
	
	SetDatafolder root:fitting:KK
	Duplicate/O w, root:fitting:KK:impart_cut
	Duplicate/O root:fitting:KK:impart_cut, root:fitting:KK:mirrorimpart_cut
	wave mirrorimpart_cut, impart_cut
	mirrorimpart_cut = 0

	variable n,nf = dimsize(impart_cut,0)
	
	for (n=0; n<(nf+1); n++)
		mirrorimpart_cut[nf-n] = impart_cut[n]
	endfor
	concatenate/NP=0/O {impart_cut, mirrorimpart_cut}, wholeIM_cut
END


// Adds a gaussian drop off for the imaginary part of the self-energy
Function gaussiancutoff(w)
	Wave w
	
	variable cutting = scaletoindex(w,0,0)	
	SetDatafolder root:fitting:KK
	
	Duplicate/O/R=[,cutting] w, w_cut1
		
	Variable startgraph = w_cut1[0]
	Variable startgraphx = 2*dimoffset(w_cut1,0)
	
	Make/O/N=(dimsize(w_cut1,0)) extension
	Setscale/P x, startgraphx, dimdelta(w_cut1,0), extension
	extension = startgraph
	
	concatenate/NP=0/O {extension, w_cut1}, w_cut

	//Variable start= -.8, stop = -0.45, width = 0.08	
	// 1 == const
	
	Variable start = -0.8, stop = -0.45, width = 0.05
	Variable middle = start-((start-stop)/2)
	variable height = w_cut[scaletoindex(w_cut,middle,0)]
	variable peak = 1/(width*sqrt(2 * pi))
	variable change = height/peak
	
	w_cut[,scaletoindex(w_cut,middle,0)] = change * Gauss(x, middle, width)
END


// Performs the Kramers-Kronig transformation (= Hilbert transform)
Function kktransform(IM)
	Wave IM
	Hilberttransform/DEST=KKIM IM
	SetScale/P x, dimoffset(IM,0), dimdelta(IM,0), KKIM
	
END


// Finds the resudal of the 2 versions of the real part
Function REresidu(REdisp, KKIM, cutoff)
	Wave REdisp, KKIM
	Variable cutoff
	
	Duplicate/O/R=[scaletoindex(KKIM,0,0),scaletoindex(KKIM,0,0) +dimsize(REdisp,0)] KKIM, KKIM_cut
	
	if (cutoff == 0)
		Duplicate/O REdisp, residual
		residual = 0
		residual = REdisp[p] - KKIM_cut[p]
	elseif (cutoff == 1)
		Duplicate/O REdisp, residualconst
		residualconst = 0
		residualconst =  REdisp[p] - KKIM_cut[p]
	endif
	
END


// The function where everything is combined to test the KK consistency for specific values
Function KKconsistent(vF,kF, EDM, cutoff)
	Variable vf, kf, cutoff
	Wave EDM
	
	fitloop(EDM)
	
	SetDataFolder root:fitting:result
	
	Wave parameterwave
	
	useparameterwave(vf,kf)
	
	Wave ReR
	
	SetDataFolder root:fitting:KK
	rawRedata(parameterwave,EDM)
	wave rawfitpeak
	
	bare_dispersion(EDM, vf, kf)
	wave dispersion

	
	ydifference(EDM)
	Wave REdisp
	
	Duplicate/O root:fitting:result:ImL, root:fitting:KK:impart
	Setdatafolder root:fitting:KK
	Wave impart
	gaussiancutoff(impart)
	Wave w_cut
	mirrorIM(w_cut)
	wave wholeIM_cut
	kktransform(wholeIM_cut)
	Wave KKIM
	
	REresidu(REdisp, KKIM, cutoff)
	
	//checkparameterwave(parameterwave,EDM, REdisp,vf,kf)
	//SetDataFolder root:fitting:KK
	
	print "2kf", abs(kf) + ReR[dimsize(ReR, 0)]

	Display KKIM, wholeIM_cut, root:KK:REdisp
	colourise(99,0)
	legend
	ModifyGraph zero(left)=1
	Label left "Sigma', Sigma\" (eV)"
	Label bottom "omega (eV)"
	
END

// Makes the 2D E, T range for the IM data
Function selfenergy2D()
	Setdatafolder root:s5:pictures
	wave s5_8K
	String folder = "root:s5:pictures"
	Variable numberOfWaves = CountObjects(folder, 1)
	Variable n
	
	Setdatafolder root:s5
	//Wave Tscale
	print numberofwaves
	Make/O/N=(dimsize(s5_8K,0), numberofwaves) selfenergy
	Setscale/P x, dimoffset(s5_8K,0), dimdelta(s5_8k,0), selfenergy

	//Setdatafolder root:s5:s5IM

	for (n = 0 ; n< numberofwaves;n++)
		print n
		Wave tmp = $(folder + ":" + GetIndexedObjName(folder, 1, n))
		selfenergy[][n] = tmp[p]
		//print n, numberofwaves-1 -n
	EndFor

END

// Makes necessary E scale for the fit
Function scaleofE(selfenergy)
	wave selfenergy
	variable l, l_final = dimoffset(selfenergy,0), l_start = indextoscale(selfenergy,66,0)
	print l_final
	Variable m =0

	Make/O/N=(dimsize(selfenergy,0)) Escale
	for (l =l_final; l<l_start+0.003; l = l +0.003)
		Escale[m] = l
		m =m+1
	endfor
	
END


// 2 dimensional T and E fit
Function fitit2D(tofit2D, withphonon)
	Wave tofit2D
	Variable withphonon
	
	SetDataFolder root:global
	If (withphonon == 1)
		Variable/G lambda=1.9, nu = 0.65, p = 1.7 , imp = 0.05, phlambda = 0.1, omegaD = 0.080 //omegaD=-0.065
	Else
		Variable/G lambda=1.9, nu = 0.7, p = 1.7, imp = 0.05 //, imp = 0.05, phlambda = 0.1, omegaD = 0.080 //omegaD=-0.065		
	Endif
	
	
	//SetDataFolder root:s5:s5IM
	//SetDataFolder root:s3IM
	//SetDataFolder root:s22
	//SetDataFolder root:s13:s13IM
	SetDataFolder root:w17prepared
	//SetDataFolder root:f5 
	//SetDataFolder root:w1b
	
	Wave Tscale, Escale
	
	
	Make/D/O/N=6 w2D = {lambda, nu, p, imp, phlambda, omegaD}
	Make/D/O/N=4 nophononw2D = {lambda, nu, p, imp}
	
	Make/O/T/N=(1) T_Constraint
	T_Constraint = {" K3 < 1.8"} //, "0.019<K4<0.02","K4>0","K5 < 0.1","K15 > 0"}
	///C=T_constraint
	
	if (withphonon ==1)
		FuncFitMD/H="010000" fitfunction2D, w2D, tofit2D[0,66][0,2] /X=Escale[0,66]/Y=Tscale[0,2]
	else
		FuncFitMD/H="0000" nophononfitfunction2D, nophononw2D, tofit2D[0,66][0,1] /C=T_constraint/X=Escale[0,66]/Y=Tscale[0,1]
	endif
	//print Tscale[2]// /C=T_constraint
	//print w2d[5]
	Variable T = 8
	if (withphonon == 1)
		if (abs(w2d[5]) > abs(dimoffset(tofit2D,0)))
			Make/O/N=(dimsize(Escale,0)) testlowT8
			Setscale/I x, -0.2, 0, testlowT8

			testlowT8 = (w2D[0] * (x*x + (w2D[2]*pi*T*w2D[2]*pi*T)/(11604.45)^2) ^(w2D[1]) +w2D[3]) + (w2d[4] *pi * -x*x*x / (3*w2d[5]*w2d[5]) )
			print "hi?"
		elseif (abs(w2d[5]) < abs(dimoffset(tofit2D,0)))
			Make/O/N=(dimsize(Escale,0)) testlowT8
			Setscale/I x, -0.2, 0, testlowT8
			PRINT "hLLO"
			Variable Ephonon = scaletoindex(tofit2D,-abs(w2d[5]),0)
			print scaletoindex(tofit2D, -0.10,0)
			testlowT8[,Ephonon] = (w2D[0] * (x*x + (w2D[2]^2*pi*T*pi*T)/(11604.45)^2) ^(w2D[1]) +w2D[3]) + (w2d[4]* w2d[5] * pi / 3) 
			testlowT8[Ephonon,] = (w2D[0] * (x*x + (w2D[2]^2*pi*T*pi*T)/(11604.45)^2) ^(w2D[1]) +w2D[3]) + (w2d[4] *pi * -x*x*x / (3*w2d[5]*w2d[5]) )

		endif
	else
		Make/O/N=(dimsize(Escale,0)) nophononlowT8
		Setscale/I x, -0.2, 0, nophononlowT8
		nophononlowT8 = (nophononw2D[0] * (x*x + (nophononw2D[2]*pi*T*nophononw2D[2]*pi*T)/(11604.45)^2) ^(nophononw2D[1]) +nophononw2D[3])
		print "hoi"
	endif
End


// fit formula for IM part
Function fitfunction2D(w2D,x,y) : FitFunc
	Wave w2D
	Variable x,y
	
	Variable lambda = w2D[0]
	Variable nu = w2D[1]
	Variable p = w2D[2]
	Variable imp = w2D[3]	
	//Variable C = w2d[4]

	
	Variable phlambda = w2D[4]
	
	Variable omegaD = w2D[5]
			
	If(abs(x) < abs(w2D[5]))
      Variable/c func = w2D[0] * (x*x + (w2D[2]*w2D[2]*pi^2*y*y)/(11604.45)^2)^(w2D[1]) + w2D[3] +( w2D[4] *pi * -x*x*x) / (3*w2D[5]*w2D[5])
 	elseif(abs(x) > abs(w2D[5]))
      func = w2D[0] * (x*x + (w2D[2]* w2D[2]*pi^2*y*y)/(11604.45)^2)^(w2D[1]) + w2D[3] + (pi * w2d[4] * w2d[5] / 3)
	endif
	
	
	return func

End


// Fit test without phonon
Function nophononfitfunction2D(w2D,x,y) : FitFunc
	Wave w2D
	Variable x,y
	
	Variable lambda = w2D[0]
	Variable nu = w2D[1]
	Variable p = w2D[2]
	Variable imp = w2D[3]	

	
	Variable func = w2D[0] * (x*x + (w2D[2]* w2D[2]*pi^2*y*y)/(11604.45)^2)^(w2D[1]) + w2D[3]
	
	return func

End


// Test phonon shape
Function testphonon()
	Variable phlambda = 0.65, omegaD = 65.
	Variable hbar = 6.582119569* 10^(-16)
	
	Make/O/N=(60) phonon
	Setscale/I x, -200, 0, phonon
	phonon =  phlambda *pi * abs(x*x*x) / (3*omegaD*omegaD)
END


// Test shape of 2D fit
Function test2Dfit(EDM)
	Wave EDM
	SetDataFolder root:s3IM
	Variable nu = 0.65, imp = 44, lambda = 0.19, p = 1.03
	Duplicate/O EDM, testselfenergy
	testselfenergy = 0
	Variable n,m
	Wave Tscale, Escale1000
	
	//print lambda * (Escale[m]*Escale[m] + (p*pi*Tscale[n]*p*pi*Tscale[n]))^(nu) + imp
	for (n = 0; n<dimsize(Tscale,0); n++)
		for (m=0; m<dimsize(Escale1000,0);m++)
			testselfenergy[m][n] = lambda * (Escale1000[m]*Escale1000[m] + (p*pi*Tscale[n]*p*pi*Tscale[n]))^(nu) + imp
		//print Tscale[n]
		endfor
	endfor
	
	Make/O/N=(dimsize(Escale1000,0)) test1D
	Setscale/I x, -200, 0, test1D
	test1D = lambda * (x*x + p*pi*Tscale[0]*p*pi*Tscale[0] ) ^(nu)	+ imp
End



//////////// Preparing EDMs //////////////


// Sets the fermi energy
Function setEf(w,Ef)
	Wave w
	Variable Ef
	
	Setscale/P y, dimoffset(w,1)-Ef, dimdelta(w,1), w
End

// k shifting and angle to k
Function shift(w,E,Enu)
	Wave w
	Variable E,Enu
	Variable xshift = getrightMDC(w,E)
	Setscale/P x, angle_to_k(dimoffset(w,0) - xshift,Enu), angle_to_k(dimdelta(w,0),Enu), w
End


Function angle_to_k(angle, Ep)
	Variable angle, Ep
	return 0.5123 * sqrt(Ep) * sin(degrees_to_rad(angle))
	
End


Function prepareEDM(w,Ef,Enu)
	Wave w
	Variable Ef,Enu
	
	shift(w,Ef,Enu)
	setEf(w,Ef)
End


Function prepareall()
	Variable/G Ef=17.6245, Enu=22
	//String folder = "root:RawData:F6ANm12"
	//String folder2 = "root:Uitwerkingen"
	String folder = "root:Fitting:Preparation:TOFIT"
	String folder2 = "root:Fitting:Prepared"
	SetDatafolder $(folder2)
	
	Variable numberOfWaves = CountObjects(folder, 1)
	Variable n
	for (n = 0; n<numberOfWaves;n++)
		Wave tmp = $(folder + ":" + GetIndexedObjName(folder, 1, n))
		Duplicate/O tmp, $(folder2 + ":" + GetIndexedObjName(folder, 1, n))
		prepareEDM($(folder2 + ":" + GetIndexedObjName(folder, 1, n)),Ef,Enu)
		print nameOfWave(tmp)
	EndFor
End


Function getEDC(w,kf)
	Wave w
	Variable kf
	
	Make/N=(dimsize(w,1)) EDC
	SetScale/P x, dimoffset(w,1), dimdelta(w,1), EDC
	EDC = w[kf][p]
End


// To correct for ARPES only accessing occupied part
Function FermiDirac(EDM, T)
	Wave EDM
	Variable T
	Variable boltzmann = 8.617 * 10^(-5)
	
	Make/O/N=(dimsize(EDM,1)) fdWave
	SetScale/P x, Dimoffset(EDM,1), Dimdelta(EDM,1), "eV", fdWave
	
	fdWave = 1 / (exp((x)/(boltzmann * T)) + 1)
End


Function dividefermi(EDM, T)
	Wave EDM
	Variable T
	
	Fermidirac(EDM, T)
	Wave fdWave
	EDM = EDM/fdWave[q]
	redimension/N=(-1,scaletoindex(EDM,0,1)) EDM
END


//F5_295K_GX = root:Fitting:F5_295K_GX[p][q] - (root:Fitting:F5_295K_GX_Hprofile0[p]/25)
//FermiDirac(F6ANm12_5K_BG)
//F6ANm12_5K_BG = F6ANm12_5K_BG/ fdWave[q]
//redimension/N=(-1,scaletoindex(F6ANm12_5K_BG,0,1)) F6ANm12_5K_BG
