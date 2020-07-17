#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


// finds an MDC at certain energy (input EDM)
Function getMDC(w,E)
	Wave w
	Variable E
	Variable Eindex = scaletoindex(w,E,1)
	
	
	Make/N=(dimsize(w,0))/O MDC
	SetScale/P x, dimoffset(w,0), dimdelta(w,0), MDC
	MDC = w[p][Eindex]
	
End

// Used to fit a lorentian over an MDC
Function fitit(tofit)
	Wave tofit
	
	SetDataFolder root:global
	NVAR a, b, c, d, e, f, g, h, i	
	SetDataFolder root:fitting:result
	
	//Make/O/N=(dimsize(tofit,0)) vorm
	//SetScale/P x, dimoffset(tofit,0), dimdelta(tofit,0), vorm
	//vorm = (1/(2 *pi)) * (a * (b/2) /((x-c)*(x-c) + (b/2)*(b/2))) +(1/(2 *pi)) * (d * (e/2) /((x-f)*(x-f) + (e/2)*(e/2) )) +g + x*h + x*x*i
	
	Make/O/N=9 w = {a, b, c, d, e, f, g, h, i}
	
	FuncFit/Q/H="000000111" fitfunction, w, tofit
End

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



Function IMfitfunction(wIM,x) : FitFunc
	Wave wIM
	Variable x 
	
	Variable T = 180
	Variable lambda = wIM[0]
	Variable nu = wIM[1]
	Variable p = wIM[2]
	Variable imp = wIM[3]	
	Variable C = wIM[4]
			
	//Impurity scattering
	Variable func = (wIM[0] * (x*x + (wIM[2]*pi*T*wIM[2]*pi*T))^(wIM[1]) + wIM[3])*wIM[4]
	
	return func
	
	
End

Function fititIM(tofit)
	Wave tofit
	
	SetDataFolder root:global
	Variable/G lambda=0.6, nu = 0.8, p = 0.9, imp = 5, C = 1
	//Variable/G nu = 0.65, imp = 5, lambda = 0.8, p = 1.03, C=1

	SetDataFolder root:s5:s5IM

	//Wave Tscale, Escale
	Wave Tscale, Escale1000
	
	Make/D/O/N=5 wIM = {lambda, nu, p, imp, C}
	
	FuncFit/Q/H="00001" IMfitfunction, wIM, tofit
	print Tscale[2]
	
	Variable T =8
	Make/O/N=(dimsize(Escale1000,0)) testfitparameters
	Setscale/I x, -0.2, 0, testfitparameters
	testfitparameters = (wIM[0] * (x*x + wIM[2]*pi*T*wIM[2]*pi*T) ^(wIM[1]) +wIM[3])*wIM[4]
	colourise(99,0)
	legend
end

//fittest = (1/ (2 *pi)) * (w[0] * w[1] / ( (x-w[2])*(x-w[2]) + (w[1]/2)*(w[1]/2) )) + (1/ (2 *pi)) * (w[3] * w[4] / ((x-w[5])*(x-w[5]) + (w[4]/2)*(w[4]/2) )) + w[6] + x * w[7] + x*x*w[8]


Function fitit2D(tofit2D, withphonon)
	Wave tofit2D
	Variable withphonon
	
	SetDataFolder root:global
	//Variable/G lambda=0.6, nu = 0.9, p = 1., imp = 0.6 , C = 1
	//Variable/G nu = 0.65, imp = 5, lambda = 0.8, p = 1.03, C=1
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
	//Wave Tscale, Escale1000
	
	//Make/D/O/N=5 w2D = {lambda, nu, p, imp, C}
	//Make/D/O/N=7 w2D = {lambda, nu, p, imp, C, phlambda, omegaD}
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
			
	//variable hbar = 1
	//Impurity scattering
	//Variable E_PH = -0.03
	
	If(abs(x) < abs(w2D[5]))
      Variable/c func = w2D[0] * (x*x + (w2D[2]*w2D[2]*pi^2*y*y)/(11604.45)^2)^(w2D[1]) + w2D[3] +( w2D[4] *pi * -x*x*x) / (3*w2D[5]*w2D[5])
 	elseif(abs(x) > abs(w2D[5]))
      func = w2D[0] * (x*x + (w2D[2]* w2D[2]*pi^2*y*y)/(11604.45)^2)^(w2D[1]) + w2D[3] + (pi * w2d[4] * w2d[5] / 3)
	endif
	//Variable func = (w2D[0] * (x*x + (w2D[2]*w2D[2]*pi^2*y*y))^(w2D[1]) + w2D[3])* w2D[4]
	return func
	
	
	
	//Sigma=cmplx(0,lam1*sqrt(E*E+(pt*pi*Tcurr)*(pt*pi*Tcurr))^(2*nu))-SE
	
	//Sigma_im = Imag(Sigma)
	return func
	
	//variable/c func =cmplx(0, (lambda*sqrt(x*x+(p*pi*y)*(p*pi*y)/(11604.45)^2)^(2*nu)+imp)*C)
	//variable Sigma_im = Imag(func)
	//return sigma_im
	
	//func += hbar * phlamba *pi * abs(x*x*x) / (3*omegaD*omegaD)
	//func = func + ( w2D[4] *pi * -x*x*x) / (3*w2D[5]*w2D[5])

End

Function nophononfitfunction2D(w2D,x,y) : FitFunc
	Wave w2D
	Variable x,y
	
	Variable lambda = w2D[0]
	Variable nu = w2D[1]
	Variable p = w2D[2]
	Variable imp = w2D[3]	

	
	Variable func = w2D[0] * (x*x + (w2D[2]* w2D[2]*pi^2*y*y)/(11604.45)^2)^(w2D[1]) + w2D[3]
	
	return func
	
	
	
	//Sigma=cmplx(0,lam1*sqrt(E*E+(pt*pi*Tcurr)*(pt*pi*Tcurr))^(2*nu))-SE
	
	//Sigma_im = Imag(Sigma)
	//return func
	
	//variable/c func =cmplx(0, (lambda*sqrt(x*x+(p*pi*y)*(p*pi*y)/(11604.45)^2)^(2*nu)+imp)*C)
	//variable Sigma_im = Imag(func)
	//return sigma_im
	
	//func += hbar * phlamba *pi * abs(x*x*x) / (3*omegaD*omegaD)
	//func = func + ( w2D[4] *pi * -x*x*x) / (3*w2D[5]*w2D[5])

End

Function fitfunctionDebye(w2D,x,y) : FitFunc
	Wave w2D
	Variable x,y
	Variable lambda = w2D[0]
	Variable boltzmann=8.617 * 10^(-5)

	//Variable func = lambda * (x) * (pi/2) * max(abs(x), boltzmann*y/(11604.45)) * sign(x)
	Variable func = lambda * abs(x)
	return func

End

//Variable Ephonon = scaletoindex(tofit2D,w2d[5],0)
	//print "test",scaletoindex(tofit2D, -0.3, 0)
	//print "Ephonon", Ephonon
	//Make/O/N=(dimsize(Escale,0)) fitT8
	//Setscale/I x, -0.2, 0, fitT8
	//Variable T = 8
	//fitT8 = (w2D[0] * (x*x + w2D[2]*pi*T*w2D[2]*pi*T/(11604.45)^2)^(w2D[1]) +w2D[3])*w2D[4]
	

	//Make/D/O/N=1 wDebye = {lambdaDebye}

	//FuncFitMD/H="1" fitfunctionDebye, wDebye, tofit2D[0,66][0,2] /X=Escale[0,66]/Y=Tscale[0,2]
	
	//Variable boltzmann=8.617 * 10^(-5) 
	//Make/O/N=(dimsize(Escale,0)) testDebye
	//Setscale/I x, -0.2, 0, testDebye
	//Variable TD = 8
	//testDebye = lambda * (x) * (pi/2) * max(abs(x), boltzmann*TD/(11604.45)) * sign(x)
	

Function testphonon()
	Variable phlambda = 0.65, omegaD = 65.
	Variable hbar = 6.582119569* 10^(-16)
	
	Make/O/N=(60) phonon
	Setscale/I x, -200, 0, phonon
	phonon =  phlambda *pi * abs(x*x*x) / (3*omegaD*omegaD)
	
END

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
ENd


// alt e o	
//	SetDataFolder root:s3IM:s3IM_cut:Fitlines
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT8
//	Setscale/I x, -200., 0, fitT8
//	fitT8 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[0]*w2D[2]*pi*Tscale[0]) ^(w2D[1]) +w2D[3])*w2D[4]
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT30
//	Setscale/I x, -200., 0, fitT30
//	fitT30 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[1]*w2D[2]*pi*Tscale[1]) ^(w2D[1]) +w2D[3])*w2D[4]
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT75
//	Setscale/I x, -200., 0, fitT75
//	fitT75 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[2]*w2D[2]*pi*Tscale[2]) ^(w2D[1]) +w2D[3])*w2D[4]
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT105
//	Setscale/I x, -200., 0, fitT105
//	fitT105 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[3]*w2D[2]*pi*Tscale[3]) ^(w2D[1]) +w2D[3])*w2D[4]
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT155
//	Setscale/I x, -200., 0, fitT155
//	fitT155 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[4]*w2D[2]*pi*Tscale[4]) ^(w2D[1]) +w2D[3])*w2D[4]
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT180
//	Setscale/I x, -200., 0, fitT180
//	fitT180 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[5]*w2D[2]*pi*Tscale[5]) ^(w2D[1]) +w2D[3])*w2D[4]
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT207
//	Setscale/I x, -200., 0, fitT207
//	fitT207 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[6]*w2D[2]*pi*Tscale[6]) ^(w2D[1]) +w2D[3])*w2D[4]
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT255
//	Setscale/I x, -200., 0, fitT255
//	fitT255 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[7]*w2D[2]*pi*Tscale[7]) ^(w2D[1]) +w2D[3])*w2D[4]
//	
//	Make/O/N=(dimsize(Escale1000,0)) fitT295
//	Setscale/I x, -200., 0, fitT295
//	fitT295 = (w2D[0] * (x*x + w2D[2]*pi*Tscale[8]*w2D[2]*pi*Tscale[8]) ^(w2D[1]) +w2D[3])*w2D[4]
//
//	
//	SetDataFolder root:s3IM


function fitloop(EDM)

	Wave EDM
	Wave parameters
	
	//Variable Ecut, E_final = scaletoindex(EDM,-1.2,1), E_start = indextoscale(EDM, dimsize(EDM, 1)-1,1)
	
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
	

	
	//Label bottom "kx(A^-1)"
	//Label left "Intensity (a.u.)"
	//SetAxis left 0,*
	//TextBox/C/N=Title/F=0/A=MC "\Z24Symmetrized EDC at kF (" + "test" + ") - Right branch"
	//TextBox/C/N=Title/A=MT/X=0.00/Y=0.00/E
	//colourise(99,0)
	//legend
	
End



// Plots values from the 2D parameterwave
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
	
	
	variable m = 208
	// Can be used to check the values of a certain slice: plots the Lorentzian
	Make/O/N=(685) check
	Setscale/I x, -0.64, 0.76, check
	check = (1/ (pi)) * (parameterwave[0][m] * (parameterwave[1][m]/2) / ( (x-parameterwave[2][m])*(x-parameterwave[2][m]) + (parameterwave[1][m]/2)*(parameterwave[1][m]/2) )) + (1/ (pi)) * (parameterwave[3][m] * (parameterwave[4][m]/2) / ((x-parameterwave[5][m])*(x-parameterwave[5][m]) + (parameterwave[4][m]/2)*(parameterwave[4][m]/2) )) + parameterwave[6][m] + x * parameterwave[7][m] + x*x*parameterwave[8][m]
	//check = (1/ (pi)) * (parameterwave[0][m] * (parameterwave[1][m]/2) / ( (x-parameterwave[2][m])*(x-parameterwave[2][m]) + (parameterwave[1][m]/2)*(parameterwave[1][m]/2) )) 

	
	//check = (1/ (pi)) * (parameterwave[0][m] * ((parameterwave[1][m]*vf)/2) / ( (x-parameterwave[2][m])*(x-parameterwave[2][m]) + ((parameterwave[1][m]*vf)/2)*((parameterwave[1][m]*vf)/2) )) 

	
	//Display check
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
		//EDM_test[][n] = (1/ pi) * ( parameterwave[0][n] *  ((parameterwave[1][n]*vf)/2) / ( (x-indextoscale(dispersion, n, 0) - REdisp[n])*(x- indextoscale(dispersion, n, 0)-REdisp[n]) + ((parameterwave[1][n]*vf)/2)*((parameterwave[1][n]*vf)/2) )) + (1/ pi) * (parameterwave[3][n] * (parameterwave[4][n]/2) / ((x-parameterwave[5][n])*(x-parameterwave[5][n]) + (parameterwave[4][n]/2)*(parameterwave[4][n]/2) )) + parameterwave[6][n] + x * parameterwave[7][n] + x*x*parameterwave[8][n]

		
		//testing fit values
		
		//EDM_values[][n] = -(1/ pi) * ( parameterwave[0][n] *  ((parameterwave[1][n]*vf)/2) / ( (x-indextoscale(dispersion, n, 0) - REdisp[n])*(x- indextoscale(dispersion, n, 0)-REdisp[n]) + ((parameterwave[1][n]*vf)/2)*((parameterwave[1][n]*vf)/2) ))
		//EDM_values[][n] = -(1/ pi) * (  ((parameterwave[1][n])/(2*vf)) / ( (x+kf-dispersion[n] +REdisp[n]/vf)*(x+kf-dispersion[n] +REdisp[n]/vf) + ((parameterwave[1][n])/2)*((parameterwave[1][n])/2) ))
		//EDM_values[][n] = -(1/ pi) * ( parameterwave[0][n] * ((parameterwave[1][n]*vf)/2) / ( (indextoscale(EDM_split,n,1) + (vf *( x - kf)) -REdisp[n])*(indextoscale(EDM_split,n,1)+ (vf* ( x - kf)) -REdisp[n]) + ((parameterwave[1][n]*vf)/2)*((parameterwave[1][n]*vf)/2) ))
		
		EDM_values[][n] = (1/ pi) * (parameterwave[0][n] * vf * ((parameterwave[1][n]*vf)/2) / ( (indextoscale(EDM_split,n,1) - (vf *( x -kf)) - REdisp[dimsize(EDM_split,1)-1 - n])^2 + ((parameterwave[1][n]*vf)/2)^2 ) )
		//EDM_values[][n] = (1/ (pi *vf)) * ((parameterwave[1][n]/(parameterwave[0][n]*2)) / ( (indextoscale(EDM_split,n,1) - (vf *( x -kf)) - REdisp[dimsize(EDM_split,1)-1 - n])^2 + ((parameterwave[1][n]/(parameterwave[0][n]*2))^2 ) ))

		
		//EDM_values[][n] = ( (1/(pi*vf))  * ((parameterwave[1][n])/2) / ((x +indextoscale(EDM_split,n,1) + kf + REdisp[dimsize(EDM_split,1)-1 - n]/vf)^2 + ((parameterwave[1][n])/2)^2 ) )
		//print REdisp[dimsize(REdisp,0)-2 -n]
		//print REdisp[dimsize(REdisp,0)-n]
		//print REdisp[dimsize(EDM_split,1)-1 - n]
	
		//EDM_values[][n] = -(1/ pi) * ( parameterwave[0][n]* ((parameterwave[1][n]*vf)/2) / ( (x - (vf * (x - kf)) -REdisp[n])*(x- (vf * (x - kf)) -REdisp[n]) + ((parameterwave[1][n]*vf)/2)*((parameterwave[1][n]*vf)/2) ))  + (1/ pi) * (parameterwave[3][n] * (parameterwave[4][n]/2) / ((x-parameterwave[5][n])*(x-parameterwave[5][n]) + (parameterwave[4][n]/2)*(parameterwave[4][n]/2) )) 
		//EDM_values[][n] = -(1/ pi) * ( parameterwave[0][n]* ((parameterwave[1][n]*vf)/2) / ( vf*vf*(x - (vf * (x - kf)) -REdisp[n])*(x- (vf * (x - kf)) -REdisp[n]) + ((parameterwave[1][n]*vf)/2)*((parameterwave[1][n]*vf)/2) ))

		
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


Function bare_dispersion(EDM,vf, kf)
	Variable vf, kf
	wave EDM
	
	SetDataFolder root:Fitting:KK
	wave rawfitpeakreversed, rawfitpeak
	
	variable last = indextoscale(EDM,0,1)/ vf + kf
	
	variable delta = -(kf - last) / dimsize(rawfitpeak,0)
	//print "delta",delta
	

	Make/O/N=(dimsize(rawfitpeak,0)) Dispersion
	//print dimsize(rawfitpeakreversed,0),rawfitpeakreversed[0]
	//dimsize(rawfitpeak,0),
	SetScale/P x, kf, delta, Dispersion
	//dimdelta(EDM,0)/1.8
	//Setscale/P y, -1.203, dimdelta(EDM,1), dispersion
	//print "deze",dimoffset(EDM,1), dimdelta(EDM,1)
	//SetScale/I y, -1.2, 0, dispersion
	//SetScale/P y, dimoffset(EDM,1), dimdelta(EDM,1), dispersion
	//dispersion = vf * (abs(x) + abs(y) - kf)
	dispersion = vf * (x - kf)

	//print rawfitpeakreversed[1] - rawfitpeakreversed[0]
	//print dispersion[0]
	//336-174 vs 3 397
	
	//calculation of delta using the rawfitpeak??

	
END

Function rawRedata(parameterwave, EDM) //, EDM)
	Wave parameterwave
	Wave EDM
	//Wave EDM
	
	SetDataFolder root:fitting:KK
	//Duplicate/O EDM, EDM_empty

	Make/O/N=(dimsize(parameterwave,1)) rawfitpeak_help
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), rawfitpeak_help
	
	
	Variable from0 = scaletoindex(EDM,0,1) -1

	rawfitpeak_help = parameterwave[2][p] //bare parameter
	Duplicate/O/R=[,from0] rawfitpeak_help, rawfitpeak

	//ModifyGraph swapXY=1

	
	//Make/O/N=(dimsize(EDM,0)) rawfitpeak2
	//SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), rawfitpeak
	//SetScale/P x, dimoffset(EDM, 0), dimdelta(EDM,0), rawfitpeak2
	//rawfitpeak2 = parameterwave[2][p]
	
	
	//Make/O/N=(dimsize(EDM,0)) rawfitpeak3
	//SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), rawfitpeak
	//SetScale/P x, dimoffset(EDM, 0), dimdelta(EDM,0), rawfitpeak3

	//rawfitpeak3[p] = rawfitpeak2[q]
	//rawfitpeak3[q] = rawfitpeak2[p]
	
	//appendtograph/VERT

	//ModifyGraph swapXY=1
	//SetAxis left -0.6,0
End

Function different()
	Wave dispersion, rawfitpeak
	variable nbegin = dimsize(rawfitpeak,0) - 1
	variable n, n_final=395

	
	variable zeropoint, m, m_final 
	
	findvalue/V=(0)/T=(0.01) dispersion
	zeropoint = V_value
	findvalue/V=(-1.2)/T=(0.01) dispersion // zelf in te stellen
	m_final = V_value
	variable dimension = m_final -zeropoint
	Make/O/N=(dimension) dispersionvalues
	//SetScale/P x, dimoffset()

	print m_final
	for (m = zeropoint; m<m_final; m++)
		//print m
		dispersionvalues[m - zeropoint] = indextoscale(dispersion,m,0)
		//print indextoscale(dispersion,m,0)
		
	endfor 
	
	reverse rawfitpeak/D=rawfitpeakreversed

	Make/O/N=(dimsize(rawfitpeakreversed,0)) dataRE
	Setscale/P x, 0, dimdelta(root:KK:LeftEDM,1), dataRE
	dataRE = -(dispersionvalues[p] - rawfitpeakreversed[p]) * pi
	//different = indextoscale(dispersion,147,0) - rawfitpeak[300]
	

End

Function ydifference(EDM)
	Wave EDM
	SetDataFolder root:Fitting:KK
	Wave dispersion, rawfitpeak, rawfitpeakreversed
	variable m = dimsize(dispersion,0), here = 300
	//print "k", rawfitpeak[m-here], "E",indextoscale(rawfitpeak, m- here,0)
	//print "E", dispersion[here], "k",indextoscale(dispersion, here,0)
	
	//print scaletoindex(dispersion,-0.3,0)
	
	
	Variable ending =  m -scaletoindex(rawfitpeak, -0.6,0)
	//print "hier",ending
	Make/O/N=(ending) REdisp
	Setscale/P x, 0, dimdelta(EDM,1), REdisp
	variable n, fitE, fitk, dispind, dispE
	for (n= 0; n< ending; n++)
		fitE = indextoscale(rawfitpeak, m-n,0)
		fitk = rawfitpeak[m-n]
		dispind = scaletoindex(dispersion, fitK,0)
		dispE = dispersion[dispind]
		//print abs(dispE) - abs(fitE)
		REdisp[n] = (abs(dispE) - abs(fitE))
		
		
	endfor
	// x waarde vinden bij y in rawfitpeak
	//print indextoscale(rawfitpeak, 372,0)
	//print rawfitpeak[372] // k waarde
	//print scaletoindex(dispersion, rawfitpeak[372],0)
	//print dispersion[54]
	
	//print abs(dispersion[54]) - abs(indextoscale(rawfitpeak, 372,0))
	
	//print indextoscale(rawfitpeak, 372,0) -
	
END


Function mirrorIM(w)
	Wave w
	
	
	SetDatafolder root:fitting:KK
	Duplicate/O w, root:fitting:KK:impart_cut
	Duplicate/O root:fitting:KK:impart_cut, root:fitting:KK:mirrorimpart_cut
	wave mirrorimpart_cut, impart_cut
	mirrorimpart_cut = 0
	//reverse/DIM=0 impart/D=mirrorimpart

	variable n,nf = dimsize(impart_cut,0)
	
	for (n=0; n<(nf+1); n++)
		mirrorimpart_cut[nf-n] = impart_cut[n]
	endfor
	
	//Make/N=(2*dimsize(impart,0)) whole_IM
	concatenate/NP=0/O {impart_cut, mirrorimpart_cut}, wholeIM_cut

END

Function gaussiancutoff(w)
	Wave w
	
	variable cutting = scaletoindex(w,0,0)
	
	//Duplicate/O w, w_cut1
	
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

Function kktransform(IM)
	Wave IM
	Hilberttransform/DEST=KKIM IM
	SetScale/P x, dimoffset(IM,0), dimdelta(IM,0), KKIM
	
END

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



Function KKconsistent(vF,kF, EDM, cutoff)
	Variable vf, kf, cutoff
	Wave EDM
	
	//fitloop(EDM)
	
	SetDataFolder root:fitting:result
	
	Wave parameterwave
	
	useparameterwave(vf,kf)
	
	Wave ReR
	
	SetDataFolder root:fitting:KK
	rawRedata(parameterwave,EDM)
	wave rawfitpeak
	
	bare_dispersion(EDM, vf, kf)
	wave dispersion
	//appendtograph/W=EDM/VERT rawfitpeak
	//appendtograph/W=EDM dispersion
	//Display
	
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
	//Setscale/P y, dimoffset(Tscale,0), dimdelta(Tscale,0), selfenergy
	//Setscale/P , 0, 300, selfenergy
	//Setscale/I x, 8,295, selfenergy
	//Label left "omega (eV)"
	//Label bottom "T"
	
	//Setdatafolder root:s5:s5IM

	for (n = 0 ; n< numberofwaves;n++)
		print n
		Wave tmp = $(folder + ":" + GetIndexedObjName(folder, 1, n))
		selfenergy[][n] = tmp[p]
		//print n, numberofwaves-1 -n
	EndFor
	
	

	//Duplicate/O/R=(l_final,-0.003)[] selfenergy, selfenergy_cut
	
	//Make/O/N=(dimsize(selfenergy,0)) Escale
	//for (l =l_final; l<l_start+0.003; l = l +0.003)
		//Escale[m] = l
		//m =m+1
	//endfor
	
END

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

//Function fitvaluetophysics(vf, kf)
	//Wave parameterwave
	//Variable vf, kf
	
	//SetDataFolder root:dispersion

//END

Function FormatGraph(datasetLabel)

	String datasetLabel

	legend
	//colourise(99, 0)
	SetAxis left -0.5,*
	//SetAxis bottom -0.2,0.2
	ModifyGraph mirror=2
	TextBox/C/N=Title/F=0/A=MC "\Z24Symmetrized EDC at kF (" + datasetLabel + ") - Right branch"
	TextBox/C/N=Title/A=MT/X=0.00/Y=0.00/E
	ModifyGraph width={Aspect,1.6},height=350
	ModifyGraph tick=2,fSize(left)=18,fSize(bottom)=18
	ModifyGraph lsize=2, width=600, height=400
	Label bottom "kx(A^-1)"
	Label left "Intensity (a.u.)"
End

Function FermiDirac(EDM, T)
	Wave EDM
	Variable T
	Variable boltzmann = 8.617 * 10^(-5)
	
	Make/O/N=(dimsize(EDM,1)) fdWave
	SetScale/P x, Dimoffset(EDM,1), Dimdelta(EDM,1), "eV", fdWave
	
	fdWave = 1 / (exp((x)/(boltzmann * T)) + 1)
	
End


////////////////////////// Part to fit things over the self energy //////////////////////////


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


Function simulationfitting(source_dispersion)
	string source_dispersion
	//NewDataFolder root:fitting:simulation
	Duplicate/O root:displayWaves:SFcut, root:fitting:simulation:$(source_dispersion)
	SetDataFolder root:fitting:simulation
	
	fitloop($(source_dispersion))
	
	NewDataFolder/O/S root:Uitwerkingen:simulation:global
			
	Variable/G j = 0.1
	
	SetDataFolder root:fitting:simulation
	
	Wave parameterwave 
	Make/O/N=(dimsize(parameterwave,1)) ImL_sim
	SetScale/P x, dimoffset(parameterwave, 1), dimdelta(parameterwave,1), ImL_sim
	ImL_sim = parameterwave[1][p] / (2 *parameterwave[0][p])  // with width W_L
	//ImL_sim = (parameterwave[1][p] * 0.2) / 2 // With some constant VF
	//ImL = parameterwave[1][p]
	Display imL_sim
	
End

Function fittingMFL(tofit)
	Wave tofit
	
	SetDataFolder root:Uitwerkingen:simulation:global
	NVAR j

	SetDataFolder root:fitting:simulation
	//Make/O/N=(dimsize(tofit,0)) vorm
	//SetScale/P x, dimoffset(tofit,0), dimdelta(tofit,0), vorm
	//vorm = (1/(2 *pi)) * (a * (b/2) /((x-c)*(x-c) + (b/2)*(b/2))) +(1/(2 *pi)) * (d * (e/2) /((x-f)*(x-f) + (e/2)*(e/2) )) +g + x*h + x*x*i
	
	Make/O/N=1 wa = {j}
	
	FuncFit/Q MFLfitfunction, wa, tofit
	Variable jj = wa[0]

	Variable T = 100, boltzmann = 8.617 * 10^(-5)
	Wave iml_sim_part
	
	Make/O/N=(100) fittest
	SetScale/P x, dimoffset(root:Fitting:simulation:IML_sim_part,0), dimdelta(root:Fitting:simulation:IML_sim_part,0), Fittest
	fittest = jj * (x) * (pi/2) * max(abs(x), boltzmann*T) * sign(x)
	Display fittest, IML_sim_part
	
End


//Function MFLfitfunction(wa,x) : FitFunc
	//Wave wa
	//Variable x
	//Variable T = 100, boltzmann = 8.617 * 10^(-5)
	
	//Variable j = wa[0] //lambda = 0.2

	//Variable func =  wa[0] * (x) * (pi/2) * max(abs(x), boltzmann*T) * sign(x)
	
	
	//return func
	
//End