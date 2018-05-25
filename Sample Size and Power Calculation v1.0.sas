***************************************************************************************************
***************************************************************************************************
**												 **
**		Sample Size and Power Calcualations for Vaccine Efficacy Trials			 **
**		via the xact Conditional and Normal Approximation Approaches			 **
**												 **
***************************************************************************************************
***************************************************************************************************
**												 **
**		Version 1.0 	Last Updated: 5/25/2018						 **
**		Developed by Matt Loiacono, MS							 **
**		Contact: MAL274@pitt.edu							 **
**												 **
**												 **
**		This code was developed as a part of a MS Thesis at the University		 **
**		of Pittsburgh Graduation School of Public Health, Department of 		 **
**		Biostatistics.									 **
**												 **
**		Please reference the following paper, by Chan & Bohidar (1998),			 **
**		for any clarifications on the formulas and methodology applied 			 **
**		in this code:									 **
**												 **
**		Chan, Ivan SF, and Norman R. Bohidar. "Exact power and sample size 		 **
**		for vaccine efficacy studies." Communications in Statistics-Theory 	 	 **
**		and Methods 27.6 (1998): 1305-1322.						 **
**												 **
***************************************************************************************************
***************************************************************************************************


**************************************************************************************************
**												**
**	MACRO: Sample Size Estimate and Power Calculation via Exact Conditional Approach 	**
**												**
**************************************************************************************************

- Calculates the exact power and estimates sample size for a vaccine efficacy trial, using an Exact 
  Conditional approach based on the Exact Conditional Test (Chan & Bohidar, 1998)
- Vaccine efficacy here is defined as (1-relative risk) of participants in the test group having the 
  disease of interest, relative to control participants having the disease
- This method defines the minimum total # of cases required (T), across the test and control group, 
  in order to achieve the pre-specified power
- Such trial is designed arounding accruing >= (T) total cases

*Required inputs are as follows (with acceptable range of values):

	P1     = disease incidence in the control group 								(0,1)
	VE1    = true vaccine efficacy of test vaccine, under the alternative hypothesis 				(0,1)
	VE0    = minimum vaccine efficacy of test vaccine, under the null (NOTE: must be less than VE1)			(0,1)
	C	   = ratio of participants (test group:control group) -- Set equal to 1 if a 1:1 ratio is desired 	(0,inf)
	Power  = pre-specificed minimum power desired to achieve							(0,1)
	Alpha  = pre-specified maximim one-sided level of the test							(0,1)
	Output = name of temporary dataset to where the results will be output

When this MACRO is run, a temporary dataset of your naming will be output, along with a printout of this dataset;	

%MACRO SampleSize_ExactConditional(P1=, VE1=, VE0=, C=, Power=, Alpha=, Output= );		
OPTIONS SPOOL;
/*Calculates and assigns theta values as macro variables for the power calculation based on input VE0 VE1 */
DATA _NULL_;
	teta0=(1-&VE0)/(1+&c-&VE0);
	teta1=(1-&VE1)/(1+&c-&VE1);
	CALL SYMPUT("theta0", teta0);
	CALL SYMPUT("theta1", teta1);
RUN;

/*Loops through all possible T's 1-2000, calculates the power and level at ALL y(obs) values*/
DATA temp;
	RETAIN VE1 P1 T Yc power level  n2;
	VE1=&VE1;
	P1=&P1;
	VE0=&VE0;
	alpha=&alpha;		*Alpha is the pre-specified significance level;
	powerreq=&power;	*powerreq is the pre-specified minimum desired power;
		DO T=1 TO 2000;
			Yc=0;
			level=0;
			power=0;
			DO x=1 TO T;
					Yc=x;	
					power = PROBBNML(&theta1,T,x);	*power is the exact power of the test;
					level = PROBBNML(&theta0,T,x);	*level is the exact level of the test;
					n2=T/((&c+1-&VE1)*(&P1));
					N=ceil(n2)*2;
					OUTPUT;
			END;
			OUTPUT;
		END;
DROP x;
RUN;

/*Removes any observations were level is > pre-specified ALPHA */
DATA temp;
	SET temp;
	IF (level <= &alpha AND power >= &power);
RUN;	

/*Keeps 1st observation PER T where minimum power is reached (this is the critical value Yc) */
DATA temp;
	SET temp;
	BY t;
	IF first.t=1;
RUN;

/* Creates FLAG1 variable where FLAG=1 if T(current)=T(previous)+1 */
DATA temp;
	SET temp;
	IF T = (LAG1(T)+1) THEN flag1=1;
	ELSE flag1=0;
RUN;	

/* Creates 20 lag variables (F1-F20) based on the FLAG variable */
DATA temp;
	SET temp;
	F1 = LAG1(FLAG1); F2 = LAG2(FLAG1); F3 = LAG3(FLAG1); F4 = LAG4(FLAG1); F5 = LAG5(FLAG1);
	F6 = LAG6(FLAG1); F7 = LAG7(FLAG1); F8 = LAG8(FLAG1); F9 = LAG9(FLAG1); F10 = LAG10(FLAG1);
	F11 = LAG11(FLAG1); F12 = LAG12(FLAG1); F13 = LAG13(FLAG1); F14 = LAG14(FLAG1); F15 = LAG15(FLAG1);
	F16 = LAG16(FLAG1); F17 = LAG17(FLAG1); F18 = LAG18(FLAG1); F19 = LAG19(FLAG1); F20 = LAG20(FLAG1);
RUN;


/* Creates new FLAG2 variable where FLAG2=1 if all flags F1 to F20 = 1, else FLAG2=0 */
/* Logic: if we count back 20 from the last FLAG2=0, then we will end up at the first T where power stabilizes (ie: doesn't drop below minimum 
	required power for any further T's)and T's become consecutive, increasing one at a time without skipping
	NOTE: due to the descrete nature of the BIN distribution, it is possible that at a given T , power is >= the pre-specified minmum power, but
	at T+1, power might drop slightly below this threshold. We want to avoid a case of acruing T+1 cases and then underpowering the trial*/
DATA temp;
	SET temp;
	IF F1=1 AND F2=1 AND F3=1 AND F4=1 AND F5=1 AND F6=1 AND F7=1 
	AND F8=1 AND F9=1 AND F10=1 AND F11=1 AND F12=1 AND F13=1 AND 
	F14=1 AND F15=1 AND F16=1 AND F17=1 AND F18=1 AND F19=1 AND F20=1 THEN FLAG2=1;
	ELSE FLAG2=0;
RUN;

/* Keeps only FLAG2=0 */
DATA temp;
	SET temp;
	WHERE FLAG2=0;
	DROP F1-F20 FLAG1;
RUN;

/*Assigns the maximum value of T in the previous datastep to MACRO variable "maxT", then subtracts 20 */
*doing this ensures we find the correct minimum value of T, based on the previous logic used;
DATA _NULL_;
	SET temp;
	CALL SYMPUT("maxT", MAX(T)-20);
RUN;

/*Identifies the lowest minimum T for which power will always be >= pre-specified power, based on maxT in previous step */
DATA &output LABEL;
	SET temp;
	IF T=&maxT;
	DROP FLAG2;
	FORMAT power level alpha powerreq PERCENT10.2;
	LABEL 	P1 = "Incidence Rate in Control"
			VE1 = "True Efficacy of Test Vaccine"
			VE0 = "Efficacy Lower Bount"
			alpha = "Pre-Specificed Level"
			powerreq = "Required Power"
			level = "Exact Level"
			power = "Exact Power"
			T = "Total # of cases (Exact Conditional)"
			N = "Expected Sample Size (Exact Conditional)";
RUN;

PROC DELETE DATA=temp;
RUN;

/* Prints output and labels accordingly */
PROC PRINT DATA=&output LABEL;
	VAR P1 VE1 VE0 powerreq alpha power level N;
RUN;

%MEND SampleSize_ExactConditional;

*Test MACRO execution for Exact Conditional Approach;
%SampleSize_ExactConditional(P1=0.05, VE1=0.65, VE0=0.15, C=1, Power=0.85, Alpha=0.05, Output=test);



***************************************************************************************************
**												 **
**	MACRO: Vaccine Efficacy Trial Sample Size using the Z Test Normal Approximation		 **
**												 **
***************************************************************************************************
-Calculates an estimated asymptotic sample size for a vaccine efficacy trial using the Noormal Approximation approach based on the Z-Test

*Required inputs are as follows (only values between 0 to 1 are accepted):


	P1     = disease incidence in the control group 								(0,1)
	VE1    = true vaccine efficacy of test vaccine, under the alternative hypothesis 				(0,1)
	VE0    = minimum vaccine efficacy of test vaccine, under the null (NOTE: must be less than VE1)			(0,1)
	C	   = ratio of participants (test group:control group) -- Set equal to 1 if a 1:1 ratio is desired 	(0,inf)
	Power  = pre-specificed minimum power desired to achieve							(0,1)
	Alpha  = pre-specified maximim one-sided level of the test							(0,1)
	Output = name of temporary dataset to where the results will be output

When this MACRO is run, a temporary dataset will be output, along with a printout of thhis dataset;	

%MACRO SampleSize_NormalApprox(P1=, VE1=, VE0=,C=, Power=, Alpha=, Output= );
OPTIONS SPOOL;

/*Calculates and assigns p2 (incidence in test group) */
DATA _NULL_;
	p2=(1-&VE1)*&P1;
	CALL SYMPUT("p2", p2);
RUN;

/* Calculates all the necessary parameters (sigma0 and sigma1) for the Normal Apprxoimation method and assigns them to macro variables */
DATA _NULL_;
	sigma1=(&P2*(1-&P2)+(1/&c)*((1-&VE0)**2)*(&P1)*(1-&P1))**(1/2);
	a=(1-&VE0)*(1+&c*&P1)+&c+&p2;
	b=(1-&VE0)*(&c*&P1+&p2);
	P2z=((a-((a**2)-4*b*(1+&c))**(1/2))/(2*(1+&c)));
	p1z=p2z/(1-&VE0);
	sigma0=(p2z*(1-p2z)+(1/&c)*((1-&VE0)**2)*(p1z)*(1-p1z))**(1/2);
	CALL SYMPUT("sigma0", sigma0);
	CALL SYMPUT("sigma1", sigma1);
RUN;

/* Calculates the estimated Normal Approximation Sample Size */
%MACRO SampleSize_NormalApprox(P1=, VE1=, VE0=,C=, Power=, Alpha=, Output= );
OPTIONS SPOOL;

/*Calculates and assigns p2 (incidence in test group) */
DATA _NULL_;
	p2=(1-&VE1)*&P1;
	CALL SYMPUT("p2", p2);
RUN;

/* Calculates all the necessary parameters (sigma0 and sigma1) for the Normal Apprxoimation method and assigns them to macro variables */
DATA _NULL_;
	sigma1=(&P2*(1-&P2)+(1/&c)*((1-&VE0)**2)*(&P1)*(1-&P1))**(1/2);
	a=(1-&VE0)*(1+&c*&P1)+&c+&p2;
	b=(1-&VE0)*(&c*&P1+&p2);
	P2z=((a-((a**2)-4*b*(1+&c))**(1/2))/(2*(1+&c)));
	p1z=p2z/(1-&VE0);
	sigma0=(p2z*(1-p2z)+(1/&c)*((1-&VE0)**2)*(p1z)*(1-p1z))**(1/2);
	CALL SYMPUT("sigma0", sigma0);
	CALL SYMPUT("sigma1", sigma1);
RUN;

/* Calculates the estimated Normal Approximation Sample Size */
DATA &output LABEL;
	P1=&p1;
	VE1=&VE1;
	VE0=&VE0;
	powerreq=&power;
	alpha=&alpha;
	zalpha=QUANTILE("norm", &alpha);
	zbeta=QUANTILE("norm", 1-&power);
	N2=ceil(((zalpha*&sigma0+zbeta*&sigma1)**2)/((&P1*(&VE1-&VE0))**2)); *sample size for test group;
	N1=ceil(N2*&c);							     *sample size for the control group;
	Ntotal=N1+N2;							     *total sample size;
	DROP zalpha zbeta;
	FORMAT power level alpha powerreq PERCENT10.2;
	LABEL 	P1 = "Incidence Rate in Control"
			VE1 = "True Efficacy of Test Vaccine"
			VE0 = "Efficacy Lower Bount"
			alpha = "Pre-Specificed Level"
			powerreq = "Required Power"
			N1 = "Excpect Sample Size in Control Group"
			N2 = "Expected Sample Size in Test Group"
			Ntotal = "Expected Total Sample Size (Normal Approximation)";
RUN;	

/* Prints output and labels accordingly */
PROC PRINT DATA=&output LABEL;
	VAR P1 VE1 VE0 powerreq alpha N1 N2 Ntotal;
RUN;

%MEND SampleSize_NormalApprox;

%MEND SampleSize_NormalApprox;

*Test MACRO execution for Normal Approximation Approach;
%SampleSize_NormalApprox(P1=0.05, VE1=0.65, VE0=0.15, C=1, Power=0.85, Alpha=0.05, Output=test);
