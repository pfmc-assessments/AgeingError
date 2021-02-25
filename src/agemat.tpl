GLOBALS_SECTION
  #include <admodel.h>
  #include <time.h>
  #include <fvar.hpp>
  dvar_vector spline(const dvector &_x,const dvar_vector&_y,dvariable yp1,
    dvariable ypn);
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
  ofstream EchoFile;

// =====================================================================

TOP_OF_MAIN_SECTION
  arrmblsize = 500000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(20000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(15000000);
  gradient_structure::set_MAX_NVAR_OFFSET(2000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(2000);
  time(&start);
  
  EchoFile.open("Echo.Out");

// =====================================================================

DATA_SECTION
 int Ireader;
 int IDataSet;
 int IDataS;
 int IDcnt;
 int MaxCells

 !! cout << "Age-reading estimation; August 2020" << endl;
 
 init_int MaxReader;
 !! EchoFile << "# Maximum readers: " << MaxReader << endl;
 init_int NDataSet;
 !! EchoFile << "# Number of data sets: " << NDataSet << endl;
 init_ivector Npnt(1,NDataSet);
 vector TotalN(1,NDataSet);
 vector EffN(1,NDataSet);
 !! EchoFile << "# Number of points per data file: " << Npnt << endl;
 init_ivector Nread(1,NDataSet);
 !! EchoFile << "# Number of readers per data file: " << Nread << endl;
 init_imatrix ReadPnt(1,NDataSet,1,Nread)
 !! EchoFile << "# Readers per data set: " << endl << ReadPnt << endl;
 init_int MinAge;
 !! EchoFile << "# Lowest age: " << MinAge << endl;
 init_int MaxAge;
 !! EchoFile << "# Maxage age: " << MaxAge << endl;
 init_ivector RefAge(1,NDataSet);
 !! EchoFile << "# RefAge (by data set): " << RefAge << endl;
 init_ivector MinusA(1,NDataSet);
 !! EchoFile << "# First age with estimated proportion (by data set): " << MinusA << endl;
 init_ivector PlusA(1,NDataSet);
 !! EchoFile << "# Last age with estimated proportion (by data set): " << PlusA << endl;
 init_ivector BiasOpt(1,MaxReader);
 !! EchoFile << "# Bias options (by READER): " << BiasOpt << endl;
 init_ivector SigOpt(1,MaxReader);
 !! EchoFile << "# Sigma options (by READER): " << SigOpt << endl;
 init_vector EffNOpt(1,NDataSet)
 !! EchoFile << "# Effective sample size by data set (<0 ignore): " << EffNOpt << endl;
 init_int UseParFile;
 !! EchoFile << "# Use PIN file (1 = Yes): " << UseParFile << endl;
 
 ivector Nreads(1,NDataSet)
 int NDataSetWithoutPerfect;                               // Number of data sets for which the answer is not known
 ivector Iperfect(1,NDataSet);                             // one if there is no perfect reading in this dataset
 !! NDataSetWithoutPerfect = 0;
 !! for (IDataS=1;IDataS<=NDataSet;IDataS++)
 !!  {
 !!   Nreads(IDataS) = 0; Iperfect(IDataS) = 0;
 !!   for (int II=1; II<=Nread(IDataS);II++)
 !!    {
 !!     if (ReadPnt(IDataS,II) > 0)
 !!      {
 !!        Nreads(IDataS) += 1;
 !!        if (SigOpt(ReadPnt(IDataS,II)) == 4) Iperfect(IDataS) = II;
 !!       }
 !!    }
 !!   if (Iperfect(IDataS) == 0) NDataSetWithoutPerfect += 1; else cout << "Age reader " << Iperfect(IDataS) << " is known age" << endl;
 !!
 !!  }
 

 // Set up bias estimation
 int NumBias;
 !! NumBias = 0; 
 !! for (Ireader=1;Ireader<=MaxReader;Ireader++)
 !! {
 !!  if (BiasOpt(Ireader) == 1) NumBias += 1;
 !!  if (BiasOpt(Ireader) == 2) NumBias += 3;
 !! }
 !! EchoFile << "Number of bias parameters " << NumBias << endl;

 // Parameters defining Bias
 vector Bias_LO(1,NumBias); 
 vector Bias_HI(1,NumBias); 
 vector Bias_INIT(1,NumBias); 
 ivector Bias_Phase(1,NumBias); 
 !! for (IDcnt=1;IDcnt<= NumBias;IDcnt++) 
 !!  *(ad_comm::global_datafile) >> Bias_LO(IDcnt) >> Bias_HI(IDcnt) >> Bias_INIT(IDcnt) >> Bias_Phase(IDcnt);
 
 // Set up Sigma estimation
 int NumSig;
 !! NumSig = 0; 
 !! for (Ireader=1;Ireader<=MaxReader;Ireader++)
 !!  {
 !!   if (SigOpt(Ireader) == 1) NumSig += 1;
 !!   if (SigOpt(Ireader) == 2) NumSig += 3;
 !!   if (SigOpt(Ireader) == 3) NumSig += 3; 
 !!  }
 
 // parameters defining Sigmas (Spline)
 matrix xvals(1,MaxReader,1,100);
 ivector nknots(1,MaxReader);
 !! for (Ireader = 1; Ireader <=MaxReader; Ireader++)
 !!  if (SigOpt(Ireader) == 5)
 !!   {
 !!    *(ad_comm::global_datafile) >> nknots(Ireader);
 !!    NumSig = NumSig + nknots(Ireader) + 2;
 !!    for (IDcnt=1;IDcnt<=nknots(Ireader);IDcnt++) *(ad_comm::global_datafile) >> xvals(Ireader,IDcnt);
 !!    EchoFile << "Spline used to define SD for reader " << Ireader << endl;
 !!    EchoFile << "Selected knots are located at ";
 !!    for (IDcnt=1;IDcnt<=nknots(Ireader);IDcnt++) EchoFile << xvals(Ireader,IDcnt) << " ";
 !!    EchoFile << endl;
 !!   } 

 // parameters defining Sigmas (Linear Interpolation)
 imatrix xvalsL(1,MaxReader,1,100);
 ivector nknotsL(1,MaxReader);
 !! for (Ireader = 1; Ireader <=MaxReader; Ireader++)
 !!  if (SigOpt(Ireader) == 6)
 !!   {
 !!    *(ad_comm::global_datafile) >> nknotsL(Ireader);
 !!    NumSig = NumSig + nknotsL(Ireader);
 !!    for (IDcnt=1;IDcnt<=nknotsL(Ireader);IDcnt++) *(ad_comm::global_datafile) >> xvalsL(Ireader,IDcnt);
 !!    if (xvalsL(Ireader,1) != 1) { cout << "First age must be 1" << endl; exit(1); }
 !!    if (xvalsL(Ireader,nknotsL(Ireader)) != MaxAge) { cout << "last age must be MaxAge" << endl; exit(1); }
 !!    EchoFile << "Linear interpolation used to define SD for reader " << Ireader << endl;
 !!    EchoFile << "Selected knots are located at ";
 !!    for (IDcnt=1;IDcnt<=nknotsL(Ireader);IDcnt++) EchoFile << xvalsL(Ireader,IDcnt) << " ";
 !!    EchoFile << endl;
 !!   } 

 !! EchoFile << "Number of variance parameters " << NumSig << endl;

 // Parameters defining Sigmas
 vector Sigma_LO(1,NumSig); 
 vector Sigma_HI(1,NumSig); 
 vector Sigma_INIT(1,NumSig); 
 ivector Sigma_Phase(1,NumSig); 
 !! for (IDcnt=1;IDcnt<= NumSig;IDcnt++) 
 !!  *(ad_comm::global_datafile) >> Sigma_LO(IDcnt) >> Sigma_HI(IDcnt) >> Sigma_INIT(IDcnt) >> Sigma_Phase(IDcnt);
    
 // Bounds and phases for the probability parameters
 init_number Prob_Low;
 init_number Prob_Hi;
 init_number Prob_Phase;

 // Determine the number of parameters that should be estimated
 int Nprobs;
 int Nslops;
 !! Nprobs = 0; Nslops = 0; 
 !! for (IDataS = 1; IDataS <=NDataSet; IDataS++)
 !!  if (Iperfect(IDataS) == 0)
 !!   {
 !!    Nprobs += PlusA(IDataS) - MinusA(IDataS);
 !!    if (MaxAge > PlusA(IDataS)) Nslops += 1;
 !!    if (MinAge < MinusA(IDataS)) Nslops += 1;
 !!   }

 // Parameters defining Slopes
 vector Slope_LO(1,Nslops); 
 vector Slope_HI(1,Nslops); 
 vector Slope_INIT(1,Nslops); 
 ivector Slope_Phase(1,Nslops); 
 !! for (IDcnt=1;IDcnt<= Nslops;IDcnt++) 
 !!  *(ad_comm::global_datafile) >> Slope_LO(IDcnt) >> Slope_HI(IDcnt) >> Slope_INIT(IDcnt) >> Slope_Phase(IDcnt);
   
 init_3darray TheData(1,NDataSet,1,Npnt,0,Nread)
 !! EchoFile << endl << "The data" << endl << TheData << endl;
 
 // -----------------------------------------------------------------------------
 
 // Find minimum and maximum ages in the data sets
 number MiAge;
 number MaAge;
 int NegVals
 !! MaAge = -1;MiAge = 1000;
 !! NegVals = 0;
 !! for (IDataS=1;IDataS<=NDataSet;IDataS++)
 !!  for (int Ipnt=1;Ipnt<=Npnt(IDataS);Ipnt++)
 !!   for (Ireader=1;Ireader<=Nread(IDataS);Ireader++)
 !!     if (TheData(IDataS,Ipnt,Ireader) >= 0)
 !!      {
 !!       if (TheData(IDataS,Ipnt,Ireader) > MaAge) MaAge = TheData(IDataS,Ipnt,Ireader);
 !!       if (TheData(IDataS,Ipnt,Ireader) < MiAge) MiAge = TheData(IDataS,Ipnt,Ireader);
 !!      }
 !!     else
 !!      NegVals = 1;
 !! if (NegVals == 1) cout << "Warning - there are some missing data; the effective sample size calculation may be dubious" << endl;

  // Storage requirements
 int NOutVal;
 !! NOutVal = (MaxReader + NumBias)*(MaxAge+1);
 
 int Nfunction;
 !! Nfunction = 0;

 ivector Ages(1,MaxReader);
 vector AgesR(1,MaxReader);
 vector x(1,10);
 int PhaseDummy;
 !! PhaseDummy = -1;

 init_int TestVal;
 !! if (TestVal != 123456) { cout << "Test Number is not 123456; it is " << TestVal << endl; exit(1); }
 !! EchoFile << TestVal << endl;

 !! EchoFile << "Structure of the data set" << endl;
 !! EchoFile << "Data set # Entries Reader boolean" << endl;
 imatrix ReaderStruc(1,1000,-1,MaxReader);
 !! ReaderStruc.initialize();
 ivector Presense(1,MaxReader);
 int NrowStruc;
 int Ifound; int Jfound;
 imatrix ReaderSumm(1,NDataSet,1,3);
 !! ReaderSumm.initialize();
 !! NrowStruc = 0;
 !! for (IDataS=1;IDataS<=NDataSet;IDataS++)
 !!  for (int II=1; II<=Npnt(IDataS);II++)
 !!   {
 !!    Presense.initialize();
 !!    for (Ireader=1;Ireader<=Nread(IDataS);Ireader++)
 !!     if (TheData(IDataS,II,Ireader) >= 0) Presense(Ireader) = ReadPnt(IDataS,Ireader);
 !!    Ifound = 0;
 !!    for (int JJ=1;JJ<=NrowStruc;JJ++)
 !!     if (ReaderStruc(JJ,-1) == IDataS)
 !!      {
 !!       Jfound = 1;
 !!       for (Ireader=1;Ireader<=Nread(IDataS);Ireader++)
 !!        if (Presense(Ireader) != ReaderStruc(JJ,Ireader)) Jfound = 0;
 !!       if (Jfound==1) Ifound = JJ;
 !!      }
 !!    if (Ifound==0)
 !!     {
 !!      NrowStruc += 1;
 !!      ReaderStruc(NrowStruc,-1) = IDataS;
 !!      ReaderStruc(NrowStruc,0) = 1;
 !!      ReaderSumm(IDataS,3) += 1;
 !!      for (Ireader=1;Ireader<=Nread(IDataS);Ireader++) ReaderStruc(NrowStruc,Ireader) = Presense(Ireader);
 !!     }
 !!    else
 !!     {
 !!      ReaderStruc(Ifound,0) += 1;
 !!     }
 !!   }
 !! for (int II=1;II<=NrowStruc;II++) EchoFile << ReaderStruc(II) << endl;
 int MaxReaderOpt;
 !! for (int II=1;II<=NrowStruc;II++)
 !!  {
 !!   IDataS = ReaderStruc(NrowStruc,-1);
 !!   ReaderSumm(IDataS,1) = IDataS;
 !!   MaxReaderOpt = 0;
 !!   for (Ireader=1;Ireader<=Nread(IDataS);Ireader++) if (ReaderStruc(NrowStruc,Ireader) >0) MaxReaderOpt += 1;
 !!   if (MaxReaderOpt > ReaderSumm(IDataS,2)) ReaderSumm(IDataS,2) = MaxReaderOpt;
 !!  }
 !!  
 !! EchoFile << "ReaderSumm" << endl;
 !! EchoFile << ReaderSumm << endl;
 
 !! cout << "Number of reads by data set:       " << Nreads << endl;
 !! cout << "Data sets without perfect readers: " << NDataSetWithoutPerfect << endl;
 !! cout << "Minimum and Maximum Ages:          " << MiAge << " " << MaAge << endl;
 
// =========================================================================
// =========================================================================

PARAMETER_SECTION
 init_bounded_number_vector SDPar(1,NumSig,Sigma_LO,Sigma_HI,Sigma_Phase);
 init_bounded_number_vector BiasPar(1,NumBias,Bias_LO,Bias_HI,Bias_Phase);
 init_bounded_vector Probs(1,Nprobs,Prob_Low,Prob_Hi,Prob_Phase)
 init_bounded_number_vector Slope(1,Nslops,Slope_LO,Slope_HI,Slope_Phase);
 init_number dummy(PhaseDummy);
 
 vector Obj_fun(1,NDataSet);
 objective_function_value obj_f;

 matrix Aprob(1,NDataSet,0,MaxAge)
 matrix TheBias(1,MaxReader,0,MaxAge)
 matrix TheSD(1,MaxReader,0,MaxAge)
 matrix Resu(1,NDataSet,1,Npnt)
 number Penal;
 
 3darray AgeErrOut(1,MaxReader,0,MaxAge,0,MaxAge)
 
 number PassIn;
 number PassOut;
 
 sdreport_vector AllOut(1,NOutVal);

// =========================================================================
// =========================================================================

PRELIMINARY_CALCS_SECTION
 dvector AgFreq(0,MaxAge);
 float Normal,Weight;
 int Age,Ipnt,AgeA,AgeB,AgeC,IDataSet,Jpnt,Iread,II;
 int JJ,Problem,OneProblem,Ifound,NLineOut,Ncells;

 // Set initial values to the values in the DAT file
 if (UseParFile != 1)
  {

   // Copy initial values
   for (JJ=1;JJ<=NumSig;JJ++) SDPar(JJ) = Sigma_INIT(JJ);
   for (JJ=1;JJ<=NumBias;JJ++) BiasPar(JJ) = Bias_INIT(JJ);
   for (JJ=1;JJ<=Nslops;JJ++) Slope(JJ) = Slope_INIT(JJ);  
   
   // Initialize the proportion parameters
   Jpnt = 0;
   for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
    if (Iperfect(IDataSet) == 0)
     {
      AgFreq.initialize();
      for (Ipnt=1;Ipnt<=Npnt(IDataSet);Ipnt++)
       for (Iread=1;Iread<=Nread(IDataSet);Iread++)
        if (TheData(IDataSet,Ipnt,Iread) >= 0)
         AgFreq(TheData(IDataSet,Ipnt,Iread)) += float(TheData(IDataSet,Ipnt,0));
      Normal = AgFreq(RefAge(IDataSet));
      for (Age=MinusA(IDataSet);Age<=PlusA(IDataSet);Age++)
       AgFreq(Age) = log( (AgFreq(Age)+0.1)/Normal);
   
      for (Age=MinusA(IDataSet);Age<=RefAge(IDataSet)-1;Age++) { Jpnt++; Probs(Jpnt) = AgFreq(Age); }
      for (Age=RefAge(IDataSet)+1;Age<=PlusA(IDataSet);Age++)  { Jpnt++; Probs(Jpnt) = AgFreq(Age); }
     }
  }   

 // Record total sample size and specify effective Ns
 for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
  {
   TotalN(IDataSet) = 0;
   for (Ipnt=1;Ipnt<=Npnt(IDataSet);Ipnt++) TotalN(IDataSet) += TheData(IDataSet,Ipnt,0);
  }
 for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
  if (EffNOpt(IDataSet) <= 0) EffN(IDataSet) = TotalN(IDataSet); else EffN(IDataSet) = EffNOpt(IDataSet);

 // Check for duplicates and condense as needed
 OneProblem = 0;
 for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
  {
   Problem = 0;
   for (II=1;II<=Npnt(IDataSet);II++)
    {
     for (JJ=1;JJ<=(II-1);JJ++)
      if (TheData(IDataSet,JJ,0) > 0)
       {
        Ifound = 0;
        for (Iread=1;Iread<=Nread(IDataSet);Iread++) if (TheData(IDataSet,JJ,Iread) != TheData(IDataSet,II,Iread)) Ifound = 1;
        if (Ifound == 0)
         {
          cout << "Warning: Lines " << II << " and " << JJ << " have the same ages" << endl; 
          TheData(IDataSet,JJ,0) += TheData(IDataSet,II,0);
          TheData(IDataSet,II,0) = -1;
          Problem = 1;
         }
       }
    }
   if (Problem == 1)
    {
     cout << "Duplicate entries found for data set " << IDataSet << "; corrected data set in Echo.Out" << endl;
     EchoFile << "Duplicate entries found for data set " << IDataSet << "; corrected data set follows" << endl;
     NLineOut = 0;
     for (II=1;II<=Npnt(IDataSet);II++)
      if (TheData(IDataSet,II,0) > 0) { NLineOut++; EchoFile << TheData(IDataSet,II) << endl; } 
     EchoFile << "New lines" << NLineOut << endl;
     OneProblem = 1; 
    }
  }
 if (OneProblem == 1) exit(1); 

 // Counter for storage
 MaxCells = 1;
 for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
  {
   Ncells = 1; for (Ireader=1;Ireader<=ReaderSumm(IDataS,2);Ireader++) { Ncells *= (MaxAge+1); }
   if (Ncells > MaxCells) MaxCells = Ncells;
  }
 cout << "total cells " << MaxCells << endl; 
 
 cout << "Done Preliminary Calcs" << endl;
 
// =========================================================================
// =========================================================================

PROCEDURE_SECTION
 int Ipnt,AgeA,Age1,IDataSet,Jreader;
 dvariable Term1,Weight,Prob1,Prob2;

 Nfunction += 1;
 
 // Initialize
 Aprob.initialize();
 Penal.initialize();
 
 // Get Bias and age-reading error CVs
 Create_Bias_and_CV();

 // Get Relative numbers by age
 Determine_Relative_Prob();

 // Create the age-reading error matrices
 Create_Error_Mats();

 // Create the aging error matrices
 obj_f = 0;
 Obj_fun.initialize();
 Resu.initialize();
 
 // Loop over all data sets
 for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
  {
 
   // Loop over points within data sets
   for (Ipnt=1;Ipnt<=Npnt(IDataSet);Ipnt++)
    {
     Weight = TheData(IDataSet,Ipnt,0)*EffN(IDataSet)/TotalN(IDataSet);
   
     // Probability associated observed data (no perfect reader)
     if (Iperfect(IDataSet) == 0)
      {
       Prob1 = 0;
       for (Age1=MinAge;Age1<=MaxAge;Age1++)
        {
         // Prior probability * product over readers
         Prob2 = Aprob(IDataSet,Age1);
         for (Ireader=1;Ireader<=Nread(IDataSet);Ireader++)
          if (int(TheData(IDataSet,Ipnt,Ireader)) >= 0)
           {
            AgeA = int(TheData(IDataSet,Ipnt,Ireader));
            Jreader = ReadPnt(IDataSet,Ireader);
            Prob2 *= AgeErrOut(Jreader,Age1,AgeA);
           } 
         Prob1 += Prob2; 
        }
      }
      
     // Probability associated observed data (a perfect reader)
     if (Iperfect(IDataSet) != 0)
      {
       // True age of this ageing structure
       Age1 = int(TheData(IDataSet,Ipnt,Iperfect(IDataSet)));
    
       // Prior probability * product over readers
       Prob1 = 1;
       for (Ireader=1;Ireader<=Nread(IDataSet);Ireader++)
        if (int(TheData(IDataSet,Ipnt,Ireader)) >= 0 & Iperfect(IDataSet) != Ireader)
         {
          AgeA = int(TheData(IDataSet,Ipnt,Ireader));
          Jreader = ReadPnt(IDataSet,Ireader);
          Prob1 *= AgeErrOut(Jreader,Age1,AgeA);
         } 
      }
     Resu(IDataSet,Ipnt) = Prob1;
     Obj_fun(IDataSet) += Weight*log(Prob1+1.0e-100);
    }
   obj_f += Obj_fun(IDataSet); 
  }
 obj_f = -1*obj_f+0.0001*Penal + 0.01*sum(square(Slope));
 obj_f += dummy*dummy;

 cout << current_phase() << " " << Nfunction << " " << obj_f << " " << Obj_fun << endl;

// ===========================================================================

FUNCTION  Create_Bias_and_CV   
 int Ireader,Ioff,Ioff3,Ioff4,IestRead,II,Iknot;
 int Age1,AgeA,BiasOption,SigOption;
 dvariable Temp,Par1,Par2,Par3,Mult,ypn,yp1;
 
 // Specify the year-specific biases
 Ioff = 1; IestRead = MaxReader-1; 
 for (Ireader=1;Ireader<=MaxReader;Ireader++) 
  {
   BiasOption = BiasOpt(Ireader);
 
   if (BiasOption == 1) 
    { Mult = BiasPar(Ioff); Ioff++; }
   if (BiasOption == 2)
    {
     Par1 = BiasPar(Ioff);
     Par2 = BiasPar(Ioff+1); 
     Par3 = BiasPar(Ioff+2);
     Temp = (Par3-Par1)/(1.0-mfexp(-Par2*(float(MaxAge)-1)));
     Ioff += 3;
    }
   for (Age1=0;Age1<=MaxAge;Age1++)
    {
     if (BiasOption == 0) TheBias(Ireader,Age1) = float(Age1)+0.5;
     if (BiasOption == 1) TheBias(Ireader,Age1) = Mult*(0.5+Age1);
     if (BiasOption == 2) TheBias(Ireader,Age1) = 0.5+Par1 + Temp*(1.0-mfexp(-Par2*(float(Age1)-1)));
     if (BiasOption == 4) TheBias(Ireader,Age1) = float(Age1)+0.5;
     if (BiasOption == 5) TheBias(Ireader,Age1) = float(Age1)+0.5;
     if (BiasOption < 0) TheBias(Ireader,Age1) = TheBias(-1*BiasOption,Age1);
    }
   // Record the estimated bias value
   if (BiasOption > 0) 
    {
     IestRead += 1;
     for (Age1=0;Age1<=MaxAge;Age1++)
      AllOut(IestRead*(MaxAge+1)+Age1+1) = TheBias(Ireader,Age1);
    }
  }  
 
 // Specify the year-specific sigmas
 Ioff = 1;
 for (Ireader=1;Ireader<=MaxReader;Ireader++) 
  {
   SigOption = SigOpt(Ireader);
 
   // Use the pattern for another reader
   if (SigOption < 0)
    {
     for (Age1=0;Age1<=MaxAge;Age1++)
      TheSD(Ireader,Age1) = TheSD(-1*SigOption,Age1);
    }

   // Parameters relate to the SD (constant CV)
   if (SigOption == 1) 
    { 
     Mult = SDPar(Ioff); 
     for (Age1=0;Age1<=MaxAge;Age1++)
      {
       if (Age1 == 0) AgeA = 1; else AgeA = Age1; 
       TheSD(Ireader,Age1) = Mult*float(AgeA);
      }
     Ioff += 1;
    }
   
   // Parameters relate to the SD (Holling Type II)
   if (SigOption == 2)
    {
     Par1 = SDPar(Ioff);
     Par2 = SDPar(Ioff+1); 
     Par3 = SDPar(Ioff+2);
     Temp = (Par3-Par1)/(1.0-mfexp(-Par2*(float(MaxAge)-1)));
     for (Age1=0;Age1<=MaxAge;Age1++)
      {
       if (Age1 == 0) AgeA = 1; else AgeA = Age1; 
       TheSD(Ireader,Age1) = Par1 + Temp*(1.0-mfexp(-Par2*(float(AgeA)-1)));
      } 
     Ioff += 3;
    }
    
   // Parameters relate to the CV
   if (SigOption == 3)
    {
     Par1 = SDPar(Ioff);
     Par2 = SDPar(Ioff+1); 
     Par3 = SDPar(Ioff+2);
     Temp = (Par3-Par1)/(1.0-mfexp(-Par2*(float(MaxAge)-1)));
     for (Age1=0;Age1<=MaxAge;Age1++)
      {
       if (Age1 == 0) AgeA = 1; else AgeA = Age1; 
       TheSD(Ireader,Age1) = (Par1 + Temp*(1.0-mfexp(-Par2*(float(AgeA)-1))))*AgeA;
      } 
     Ioff += 3;
    }
   
   // No error
   if (SigOption == 4)
    for (Age1=0;Age1<=MaxAge;Age1++) TheSD(Ireader,Age1) = 0;

   // Spline Parameters relate to the CV
   if (SigOption == 5)
    {
     cout << "Spline" << endl;
     dvector x(1,nknots(Ireader));
     dvar_vector VyI(1,nknots(Ireader));
     for (II = 1; II<=nknots(Ireader);II++) 
      { x(II) = xvals(Ireader,II); VyI(II) = mfexp(SDPar(Ioff+II-1)); } 
     Ioff += nknots(Ireader);
     yp1 = SDPar(Ioff);
     ypn = SDPar(Ioff+1);
     Ioff += 2;
     dvar_vector ders = spline(x,VyI,yp1,ypn);
     for (Age1=0;Age1<=MaxAge;Age1++)
      {
       if (Age1 == 0) AgeA = 1; else AgeA = Age1; 
       TheSD(Ireader,Age1) = sqrt(square(splint(x,VyI,ders,float(AgeA))));
      } 
    }
   
   // Linear interpolation
   if (SigOption == 6)
    {
     dvector x(1,nknotsL(Ireader));
     dvar_vector VyI(1,nknotsL(Ireader));
     for (II = 1; II<=nknotsL(Ireader);II++) 
      { x(II) = float(xvalsL(Ireader,II)); VyI(II) = mfexp(SDPar(Ioff+II-1)); } 
     Ioff += nknotsL(Ireader);
     for (Iknot=1;Iknot<=nknotsL(Ireader)-1;Iknot++)
      {
       Par1 = (VyI(Iknot+1)-VyI(Iknot))/(x(Iknot+1)-x(Iknot));
       Par2 = VyI(Iknot) - Par1*x(Iknot);
       for (Age1=xvalsL(Ireader,Iknot);Age1<xvalsL(Ireader,Iknot+1);Age1++)
        TheSD(Ireader,Age1) = Par2 + Par1*float(Age1);
      }
     TheSD(Ireader,MaxAge) = Par2 + Par1*float(MaxAge); 
     TheSD(Ireader,0) = TheSD(Ireader,1);
    }
    
   // Save results
   for (Age1=0;Age1<=MaxAge;Age1++)
    AllOut((Ireader-1)*(MaxAge+1)+Age1+1) = TheSD(Ireader,Age1);
     
  }  
   
// ===========================================================================

FUNCTION cumd_norm2
 dvariable tmp;
 
// if (PassIn < -20.0)
//  PassOut = 1.0e-190;
// else
//  if (PassIn > 20.0)
//   PassOut = 1;
//  else
//   PassOut = cumd_norm(PassIn);
   tmp = -20 + 40/(1+mfexp(-0.1*PassIn));
   PassOut = cumd_norm(tmp);
 
// ===========================================================================
// ===========================================================================

FUNCTION Determine_Relative_Prob
 dvariable Term1,AllTotal,Total,Temp;
 int Age1,Ioff,Joff,AgeA,AgeB,AgeC,IDataSet;
 int Alow,Ahi,Alow1,Ahi1,Alow2,Ahi2,Alow3,Ahi3;
 int Iread1, Iread2, Iread3;

 Ioff = 1; Joff = 1;
 Aprob.initialize();
 for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
  {
 
   // Only deal with this for cases where this no perfect reader
   if (Iperfect(IDataSet) == 0)
    {
     Total = 0;
     for (Age1=MinusA(IDataSet);Age1<=PlusA(IDataSet);Age1++)
      {
       if (Age1 != RefAge(IDataSet))
        {
         Aprob(IDataSet,Age1) = mfexp(Probs(Ioff));
         Penal += Probs(Ioff)*Probs(Ioff);
         Total += Aprob(IDataSet,Age1);
         Ioff++;
        }
       else
        {
         Aprob(IDataSet,RefAge(IDataSet)) = 1.0;
         Total += 1.0;
        } 
      }
     for (Age1=MinAge;Age1<MinusA(IDataSet);Age1++)
      {
       Temp = Slope(Joff)*(MinusA(IDataSet)-Age1);
       Aprob(IDataSet,Age1) = Aprob(IDataSet,MinusA(IDataSet))*mfexp(Temp);
       Total += Aprob(IDataSet,Age1);
      }
     if (MinusA(IDataSet) > MinAge) Joff += 1;
     for (Age1=PlusA(IDataSet)+1;Age1<=MaxAge;Age1++)
      {
       Temp = Slope(Joff)*(Age1-PlusA(IDataSet));
       Aprob(IDataSet,Age1) = Aprob(IDataSet,PlusA(IDataSet))*mfexp(Temp);
       Total += Aprob(IDataSet,Age1);
      }
     if (PlusA(IDataSet) < MaxAge) Joff += 1;
      
     // Normalize 
     for (Age1=0;Age1<=MaxAge;Age1++) Aprob(IDataSet,Age1) /= Total;
    }  
 
  }
 
// ===========================================================================

FUNCTION Create_Error_Mats
  int Ireader,Age1,Age2;
  dvariable Temp, Diff,Total, SDD1,Prop,tot;
  dvar_vector tmp1(0,MaxAge);

  // Matrices are constructed to match Synthesis
  for (Ireader=1;Ireader<=MaxReader;Ireader++) 
   for (Age1=0;Age1<=MaxAge;Age1++)
    {
     tmp1.initialize(); tot=0;
     SDD1 = TheSD(Ireader,Age1)+1.0e-30;
     for (Age2=1;Age2<=MaxAge;Age2++)
      {
       Diff = float(Age2) - TheBias(Ireader,Age1);
       PassIn = Diff/SDD1; cumd_norm2(); tmp1(Age2-1) = PassOut-tot;
       tot = PassOut;
      }  
     tmp1(MaxAge) = 1.0 - tot; 
     Total = 0;
     for (Age2=0;Age2<=MaxAge;Age2++)
      {
       AgeErrOut(Ireader,Age1,Age2) = tmp1(Age2);
       Total += tmp1(Age2);
      }
    for (Age2=0;Age2<=MaxAge;Age2++) AgeErrOut(Ireader,Age1,Age2) /= Total;
   }  

// ===========================================================================

FINAL_SECTION
 time(&finish); 
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;

// ===========================================================================

REPORT_SECTION
 int Age,Ipnt,IDataSet,Ireader,AgeA,AgeB,AgeC,Iread,Jread,Incc;
 int Ncells,Ndiv,DivI,DivJ,Age1,II,JJ,KK,Ifound,Jfound,Jreader;
 dvariable CVV,EffPred,Nages,TotalProb,Prob1,Prob2,Pobs,Pest,Top,Bot;
 dvar_vector ProbStore(1,MaxCells);
 
 report << obj_f << " " << Obj_fun << endl;
 report << "Total number of readers: " << MaxReader << endl;
 report << "Number of data sets: " << NDataSet << endl;
 report << "Bias options by reader: " << BiasOpt << endl;
 report << "Sigma options by reader: " << SigOpt << endl;
 
 report << "Total objective function: " << obj_f << endl;
 report << "Number of readers: " << MaxReader << endl;
 report << "Range of ages: " << MinAge << " - " << MaxAge << endl;
 for (Ireader=1;Ireader<=NDataSet;Ireader++)
  report << "Reader #" << Ireader << " Minus/Plus ages: " << MinusA(Ireader) << " / " << PlusA(Ireader) << endl;
 report << "Number of data sets: " << NDataSet << endl;
 report << "Number of lines of data per data set: " << Npnt << endl;
 report << "Number of data points per data set:" << TotalN << endl;
 report << endl;
 report << "Reader BiasType SigmaType" << endl;
 for (Ireader=1;Ireader<=MaxReader;Ireader++)
  report << Ireader << " " << BiasOpt(Ireader) << " " << SigOpt(Ireader) << endl;
 report << endl; 

 // Bias and variance 
 report << "Reader Age CV SD Expected age" << endl;
 for (Ireader=1;Ireader<=MaxReader;Ireader++)
  for (Age=0;Age<=MaxAge;Age++)
   {
    if (Age > 1)
     CVV =  TheSD(Ireader,Age)/float(Age);
    else
     CVV = TheSD(Ireader,Age); 
   report << Ireader << " " << Age << " " << CVV << " " << TheSD(Ireader,Age)  << " " << TheBias(Ireader,Age) << endl;
  }
 report << endl; 
 
 // Estimated age-structure
 report << "Estimated age-structure by data set" << endl;
 report << "===================================" << endl;
 report << "Age "; for (IDataSet=1;IDataSet<=NDataSet;IDataSet++) report << "Data set#" << IDataSet << " "; report << endl;
 for (Age=0;Age<=MaxAge;Age++)
  {
   report << Age << " ";
   for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
    report <<  Aprob(IDataSet,Age) << " ";
   report << endl;
  }  
 report << endl; 
 
 report << "Final age-reading error matrices" << endl;
 for (Ireader=1;Ireader<=MaxReader;Ireader++) 
  {
   report << "Matrix for reader# " << Ireader << endl;
   report << AgeErrOut(Ireader) << endl;
  }
 report << endl; 
 
 // find the effective sample sizes
 report << "Compute the effective sample sizes" << endl;
 report << "==================================" << endl;
 int GroupPointer;
 GroupPointer = 0;
 if (last_phase())
 for (IDataSet=1;IDataSet<=NDataSet;IDataSet++)
  {
   report << "Data set: " << IDataSet << endl;
   report << "Data_set Group Group Line Readers Obs Obs_probability Pred_Probability" << endl;
   Top = 0; Bot = 0; 
   
   for (int Kgroup=1;Kgroup<=ReaderSumm(IDataSet,3);Kgroup++)
    {
     GroupPointer += 1;
   
     // Find the total number of combinations of ages (MaxAge+1)**number_of_readers
     Ncells = 1; for (Iread=1;Iread<=ReaderSumm(IDataS,2);Iread++) Ncells *= (MaxAge+1);
   
     // Move through each possible combination of ages between 0 and MaxAge 
     TotalProb = 0;
     for (II=1;II<=Ncells;II++)
      {
     
       // Find the ages for this "cell" 
       Ndiv = II;
       for (Iread=1;Iread<=Nread(IDataSet);Iread++) Ages(Iread) = -1;
       for (Iread=1;Iread<=Nread(IDataSet);Iread++)
        if (ReaderStruc(GroupPointer,Iread) > 0)
         {
          DivJ = 1; 
          for (Jread=Nread(IDataSet);Jread>=Iread+1;Jread--) 
           if (ReaderStruc(GroupPointer,Jread) > 0) DivJ *= (MaxAge+1);
          DivI = (Ndiv-1)/DivJ;    
          Ages(Iread) = DivI;
          Ndiv = Ndiv - DivI*DivJ;
         }
      
      // Find the probability for this cell, i.e. the probability of an ageing reading of Ages(1)&Ages(2)&...
      if (Iperfect(IDataSet) == 0)
       {
        Prob1 = 0;
        for (Age1=MinAge;Age1<=MaxAge;Age1++)
         {
          // Prior probability * product over readers
          Prob2 = Aprob(IDataSet,Age1);
          for (Ireader=1;Ireader<=Nread(IDataSet);Ireader++)
           { 
            Jreader = ReadPnt(IDataSet,Ireader);
            AgeA = Ages(Ireader); 
            if (AgeA >= 0) Prob2 *= AgeErrOut(Jreader,Age1,AgeA); 
           } 
          Prob1 += Prob2; 
         }
        ProbStore(II) = Prob1; 
        TotalProb += Prob1;   
       }
      else
       {
        // Product over readers
        Age1 = Ages(1);
        Prob1 = 1;
        for (Ireader=2;Ireader<=Nread(IDataSet);Ireader++)
         { 
          Jreader = ReadPnt(IDataSet,Ireader);
          AgeA = Ages(Ireader); 
          if (Age1 >=0 & AgeA >= 0) Prob1 *= AgeErrOut(Ireader,Age1,AgeA); 
         }
        ProbStore(II) = Prob1; 
        TotalProb += Prob1;   
       }
      }
     
     // Now complute 
     for (II=1;II<=Ncells;II++)
      {
      
       // Find the ages for this "cell" 
       Ndiv = II;
       for (Iread=1;Iread<=Nread(IDataSet);Iread++) Ages(Iread) = -1;
       for (Iread=1;Iread<=Nread(IDataSet);Iread++)
        if (ReaderStruc(GroupPointer,Iread) > 0)
         {
          DivJ = 1; 
          for (Jread=Nread(IDataSet);Jread>=Iread+1;Jread--) 
           if (ReaderStruc(GroupPointer,Jread) > 0) DivJ *= (MaxAge+1);
          DivI = (Ndiv-1)/DivJ;    
          Ages(Iread) = DivI;
          Ndiv = Ndiv - DivI*DivJ;
         }
      
       // Check for a match
       Jfound = 0; Pobs = 0;
       for (JJ=1;JJ<=Npnt(IDataSet);JJ++)
        {
         // Set Ifound to 1 if the ages don't match properly
         Ifound = 0;
         for (Iread=1;Iread<=Nread(IDataSet);Iread++)
          if (Ages(Iread) != int(TheData(IDataSet,JJ,Iread)) & int(TheData(IDataSet,JJ,Iread)) >=0) Ifound = 1;
         // We have a match so store the results 
         if (Ifound == 0)
          {  Pobs = TheData(IDataSet,JJ,0)/ReaderStruc(GroupPointer,0); Jfound = 1; KK = JJ; }
         }
      
       // Find the probability for this cell 
       Pest = ProbStore(II) / TotalProb+1.0e-190; 
       if (Jfound == 1) 
        {
         report << "Data Point: " << IDataSet << " " << GroupPointer << " " << Kgroup << " " << KK << " ";
         for (Iread=1;Iread<=Nread(IDataSet); Iread++) report << Ages(Iread) << " ";
         report << TheData(IDataSet,KK,0) << " " << Pobs << " " << Pest << " " << ReaderStruc(GroupPointer,0) << endl;
        } 
      
       // Compute the effective sample size muliplier
       Top += Pest*(1-Pest)/float(ReaderStruc(GroupPointer,0));
       Bot += square(Pest-Pobs);
      
      }
     }
    EffPred = Top/Bot*TotalN(IDataSet);
    report << "Data_Set Predicted_EFF Assumed_Eff Sample_size" << endl;
    report << IDataSet << " " << EffPred << " " << EffN(IDataSet) << " " << TotalN(IDataSet) << endl;
    if (EffNOpt(IDataSet) < 0) EffN(IDataSet) = value(EffPred);
    report << endl;
   }
  cout << "Completed the Report Section" << endl;


// ===========================================================================

RUNTIME_SECTION
  convergence_criteria 1.e-3,1.e-6,1.e-9
  maximum_function_evaluations 1000, 50000
