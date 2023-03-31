#define TMB_LIB_INIT R_init_AgeingError
#include <TMB.hpp>

template <class Type>
struct general{
  vector<int> BiasOpt;                                                        // Bias specifications by reader
  vector<int> SigOpt;                                                         // SD specifications by reader
  vector<int> RefAge;                                                         // Reference age by data set
  vector<int> MinusA;                                                         // Minus group
  vector<int> PlusA;                                                          // Plus group
  vector<int> Iperfect;                                                       // Is the current reader perfect (1) or not (0)
  matrix<int> xvals;                                                          // x-values for spline model
  vector<int> nknots;                                                         // number of knots (spline model)
  matrix<int> xvalsL;                                                         // x-values for linear model
  vector<int> nknotsL;                                                        // number of knots (linear model)
  int NDataSet;
  int MaxReader;
  int MaxAge;
  int MinAge;
};



//=============================================================================================

template <class Type>
 Type Create_Bias_and_CV(general<Type> &gen, vector<Type> BiasPar, vector<Type> SDPar, matrix<Type> &TheBias, matrix<Type> &TheSD) {
   int Ioff, BiasOption, SigOption, AgeA;
   Type Mult, Par1, Par2, Par3, Temp;

   // Specify the year-specific sigmas bias
   Ioff = 0;
   for (int Ireader=0;Ireader<gen.MaxReader;Ireader++)
	{
     BiasOption = gen.BiasOpt(Ireader);

     // Constant bias
     if (BiasOption == 1)
      { Mult = BiasPar(Ioff); Ioff++; }
     // Asymptotic bias
     if (BiasOption == 2)
      {
       Par1 = BiasPar(Ioff);
       Par2 = BiasPar(Ioff+1);
       Par3 = BiasPar(Ioff+2);
       Temp = (Par3-Par1)/(1.0-exp(-Par2*(float(gen.MaxAge)-1)));
       Ioff += 3;
      }
     // Made prediction
     for (int Age1=0;Age1<=gen.MaxAge;Age1++)
      {
       if (BiasOption == 0) TheBias(Ireader,Age1) = float(Age1)+0.5;
       if (BiasOption == 1) TheBias(Ireader,Age1) = Mult*(0.5+Age1);
       if (BiasOption == 2) TheBias(Ireader,Age1) = 0.5+Par1 + Temp*(1.0-exp(-Par2*(float(Age1)-1)));
       if (BiasOption == 4) TheBias(Ireader,Age1) = float(Age1)+0.5;
       if (BiasOption == 5) TheBias(Ireader,Age1) = float(Age1)+0.5;
       if (BiasOption < 0) TheBias(Ireader,Age1) = TheBias(-1*BiasOption-1,Age1);            // -1 because arrays start at zero
      }
	} // Ireader

   // Specify the year-specific sigmas
   Ioff = 0;
   for (int Ireader=0;Ireader<gen.MaxReader;Ireader++)
    {
     SigOption = gen.SigOpt(Ireader);

    // Use the pattern for another reader
    if (SigOption < 0)
     {
      for (int Age1=0;Age1<=gen.MaxAge;Age1++)
       TheSD(Ireader,Age1) = TheSD(-1*SigOption-1,Age1);                                    // -1 because arrays start at zerp
     }

    // Parameters relate to the SD (constant CV)
    if (SigOption == 1)
     {
      Mult = SDPar(Ioff);
      for (int Age1=0;Age1<=gen.MaxAge;Age1++)
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
      Temp = (Par3-Par1)/(1.0-exp(-Par2*(float(gen.MaxAge)-1)));
      for (int Age1=0;Age1<=gen.MaxAge;Age1++)
       {
        if (Age1 == 0) AgeA = 1; else AgeA = Age1;
        TheSD(Ireader,Age1) = Par1 + Temp*(1.0-exp(-Par2*(float(AgeA)-1.0)));
        }
      Ioff += 3;
     }

    // Parameters relate to the CV
    if (SigOption == 3)
     {
      Par1 = SDPar(Ioff);
      Par2 = SDPar(Ioff+1);
      Par3 = SDPar(Ioff+2);
      Temp = (Par3-Par1)/(1.0-exp(-Par2*(float(gen.MaxAge)-1)));
      for (int Age1=0;Age1<=gen.MaxAge;Age1++)
       {
        if (Age1 == 0) AgeA = 1; else AgeA = Age1;
        TheSD(Ireader,Age1) = (Par1 + Temp*(1.0-exp(-Par2*(float(AgeA)-1.0))))*AgeA;
       }
      Ioff += 3;
      }

    // No error
    if (SigOption == 4)
     for (int Age1=0;Age1<=gen.MaxAge;Age1++) TheSD(Ireader,Age1) = 0;

    // Spline Parameters relate to the CV
    if (SigOption == 5)
     {
      Type Temp2;
      int AgeA;
      vector<Type>  x(gen.nknots(Ireader));
      vector<Type> VyI(gen.nknots(Ireader));
      for (int II = 0; II<gen.nknots(Ireader);II++)
       { x(II) = float(gen.xvals(Ireader,II)); VyI(II) = exp(SDPar(Ioff+II)); }
      Ioff += gen.nknots(Ireader);
      tmbutils::splinefun<Type> spline(x,VyI,3);
      for (int Age1=0;Age1<=gen.MaxAge;Age1++)
       {
        if (Age1 == 0) AgeA = 1; else AgeA = Age1;
        Temp2 = float(AgeA);
        Temp2 = spline(Temp2);
        TheSD(Ireader,Age1) = sqrt(Temp2*Temp2);
       }
	 }

    // Linear interpolation
    if (SigOption == 6)
     {
      Type Par1; Type Par2;
      vector<Type> x(gen.nknotsL(Ireader));
      vector<Type> VyI(gen.nknotsL(Ireader));
      for (int II = 0; II<gen.nknotsL(Ireader);II++)
       { x(II) = float(gen.xvalsL(Ireader,II)); VyI(II) = exp(SDPar(Ioff+II)); }
      Ioff += gen.nknotsL(Ireader);
      for (int Iknot=0;Iknot<(gen.nknotsL(Ireader)-1);Iknot++)
       {
        Par1 = (VyI(Iknot+1)-VyI(Iknot))/(x(Iknot+1)-x(Iknot));
        Par2 = VyI(Iknot) - Par1*x(Iknot);
        for (int Age1=gen.xvalsL(Ireader,Iknot);Age1<gen.xvalsL(Ireader,Iknot+1);Age1++)
         TheSD(Ireader,Age1) = Par2 + Par1*float(Age1);
       }
      TheSD(Ireader,gen.MaxAge) = Par2 + Par1*float(gen.MaxAge);
      TheSD(Ireader,0) = TheSD(Ireader,1);
	 }

	// Linear change in SD with age
	if (SigOption == 7)
	 {
      Par1 = SDPar(Ioff);
      Par2 = SDPar(Ioff+1);
      Temp = (Par2-Par1)/(float(gen.MaxAge)-1);
      for (int Age1=0;Age1<=gen.MaxAge;Age1++)
       {
        if (Age1 == 0) AgeA = 1; else AgeA = Age1;
        TheSD(Ireader,Age1) = (Par1 + Temp*(float(AgeA)-1.0));
       }
      Ioff += 2;
	 }

	// Linear change in CV with age
	if (SigOption == 8)
	 {
      Par1 = SDPar(Ioff);
      Par2 = SDPar(Ioff+1);
      Temp = (Par2-Par1)/(float(gen.MaxAge)-1);
      for (int Age1=0;Age1<=gen.MaxAge;Age1++)
       {
        if (Age1 == 0) AgeA = 1; else AgeA = Age1;
        TheSD(Ireader,Age1) = (Par1 + Temp*(float(AgeA)-1.0))*AgeA;
       }
      Ioff += 2;
	 }

   } // Ireader

  Type XX;
  XX = 5.0;
  return XX;
 }

// -----------------------------------------------------------------------------

template <class Type>
 Type Determine_Relative_Prob(general<Type> &gen, vector<Type> Probs, vector<Type> Slope, matrix<Type> &Aprob, Type &Penal) {
  int Ioff, Joff;
  Type Total, Temp;

  Ioff = 0; Joff = 0;
  Aprob.setZero();
  Penal = 0;
  for (int IDataSet=0;IDataSet<gen.NDataSet;IDataSet++)
   {
   // Only deal with this for cases where this no perfect reader
   if (gen.Iperfect(IDataSet) == 0)
    {
     Total = 0;
     for (int Age1=gen.MinusA(IDataSet);Age1<=gen.PlusA(IDataSet);Age1++)
      {
       if (Age1 != gen.RefAge(IDataSet))
        {
         Aprob(IDataSet,Age1) = exp(Probs(Ioff));
         Penal += Probs(Ioff)*Probs(Ioff);
         Total += Aprob(IDataSet,Age1);
         Ioff++;
        }
       else
        {
         Aprob(IDataSet,gen.RefAge(IDataSet)) = 1.0;
         Total += 1.0;
        }
      }

     // Fill in below minus and above plus-group
     for (int Age1=gen.MinAge;Age1<gen.MinusA(IDataSet);Age1++)
      {
       Temp = Slope(Joff)*(gen.MinusA(IDataSet)-Age1);
       Aprob(IDataSet,Age1) = Aprob(IDataSet,gen.MinusA(IDataSet))*exp(Temp);
       Total += Aprob(IDataSet,Age1);
      }
     if (gen.MinusA(IDataSet) > gen.MinAge) Joff += 1;
     for (int Age1=gen.PlusA(IDataSet)+1;Age1<=gen.MaxAge;Age1++)
      {
       Temp = Slope(Joff)*(Age1-gen.PlusA(IDataSet));
       Aprob(IDataSet,Age1) = Aprob(IDataSet,gen.PlusA(IDataSet))*exp(Temp);
       Total += Aprob(IDataSet,Age1);
      }
     if (gen.PlusA(IDataSet) < gen.MaxAge) Joff += 1;

     // Normalize
     for (int Age1=0;Age1<=gen.MaxAge;Age1++) Aprob(IDataSet,Age1) /= Total;
    }
   }

  Type XX;
  XX = 5.0;
  return XX;
 }

// -----------------------------------------------------------------------------

template <class Type>
 Type Create_Error_Mats(general<Type> &gen, array<Type> &AgeErrOut,matrix<Type> &TheBias, matrix<Type> &TheSD, matrix<Type> &Aprob) {
  vector<Type> tmp1(gen.MaxAge+1);
  Type tot,SDD1,Diff,PassOut,Total,PassIn,tmp;

  // Matrices are constructed to match Synthesis
  for (int Ireader=0;Ireader<gen.MaxReader;Ireader++)
   for (int Age1=0;Age1<=gen.MaxAge;Age1++)
    {
     tmp1.setZero(); tot=0;
     SDD1 = TheSD(Ireader,Age1)+1.0e-30;
     for (int Age2=1;Age2<=gen.MaxAge;Age2++)
      {
       Diff = float(Age2) - TheBias(Ireader,Age1);
       PassIn = Diff/SDD1;
       tmp = -20.0 + 40.0/(1+exp(-0.1*PassIn));
       PassOut = pnorm(tmp,Type(0),Type(1));
       tmp1(Age2-1) = PassOut-tot;
       tot = PassOut;
      }
     tmp1(gen.MaxAge) = 1.0 - tot;
     Total = 0;
     for (int Age2=0;Age2<=gen.MaxAge;Age2++)
      {
       AgeErrOut(Ireader,Age1,Age2) = tmp1(Age2);
       Total += tmp1(Age2);
      }
    for (int Age2=0;Age2<=gen.MaxAge;Age2++) AgeErrOut(Ireader,Age1,Age2) /= Total;
   }

  Type XX;
  XX = 5.0;
  return XX;
 }

// ==============================================================================

template<class Type>
Type objective_function<Type>::operator() ()
{
 general<Type> General;
 DATA_INTEGER(NDataSet); General.NDataSet = NDataSet;
 DATA_INTEGER(MinAge); General.MinAge = MinAge;
 DATA_INTEGER(MaxAge); General.MaxAge = MaxAge;
 DATA_IVECTOR(BiasOpt); General.BiasOpt = BiasOpt;
 DATA_IVECTOR(SigOpt); General.SigOpt = SigOpt;
 DATA_INTEGER(MaxReader); General.MaxReader = MaxReader;
 DATA_IVECTOR(MinusA); General.MinusA = MinusA;
 DATA_IVECTOR(PlusA); General.PlusA = PlusA;
 DATA_IVECTOR(RefAge); General.RefAge = RefAge;
 DATA_IVECTOR(Iperfect); General.Iperfect = Iperfect;
 DATA_IVECTOR(nknots); General.nknots = nknots;
 DATA_IMATRIX(xvals); General.xvals = xvals;
 DATA_IVECTOR(nknotsL); General.nknotsL = nknotsL;
 DATA_IMATRIX(xvalsL); General.xvalsL = xvalsL;
 DATA_IARRAY(TheData);
 DATA_IVECTOR(Npnt);
 DATA_IVECTOR(Nread);
 DATA_INTEGER(MaxNpnt);
 DATA_IMATRIX(ReadPnt);
 DATA_SCALAR(AprobWght);
 DATA_SCALAR(SlopeWght);

 PARAMETER_VECTOR(BiasPar);
 PARAMETER_VECTOR(SDPar);
 PARAMETER_VECTOR(Slope);
 PARAMETER_VECTOR(Probs);
 PARAMETER(Dummy);

 matrix<Type> TheBias(MaxReader,MaxAge+1);                // Bias (note MaxAge+1 because we start at age zero)
 TheBias.setZero();
 matrix<Type> TheSD(MaxReader,MaxAge+1);                  // Standard deviation (note MaxAge+1 because we start at age zero)
 TheSD.setZero();
 matrix<Type> Aprob(NDataSet,MaxAge+1);                   // Relative probabilities
 Aprob.setZero();
 array<Type>AgeErrOut(MaxReader,MaxAge+1,MaxAge+1);       // Ageing error matrix
 AgeErrOut.setZero();

 matrix<Type> Resu(NDataSet,MaxNpnt);                     // Results matrix
 vector<Type> Obj_fun(NDataSet);                          // Objective function by data set
 Type Weight,Prob1,Prob2;                                 // Temps
 Type f;                                                  // Objective function
 Type Penal,PenalSlope;                                   // Penalty function
 Type XX;                                                 // Passing variable
 int AgeA,AgeTrue,Jreader;


 // Get the bias and age-reading error CVs
 XX = Create_Bias_and_CV(General,BiasPar,SDPar,TheBias,TheSD);

 // Determine the relative probability by age-class
 XX = Determine_Relative_Prob(General, Probs, Slope, Aprob, Penal);

 // Compute the error matrices
 XX = Create_Error_Mats(General,AgeErrOut,TheBias,TheSD,Aprob);

 // Loop over all data sets
 f = 0;
 Obj_fun.setZero();
 for (int IDataSet=0;IDataSet<NDataSet;IDataSet++)
  {

   // Loop over points within data sets
   for (int Ipnt=0;Ipnt<Npnt(IDataSet);Ipnt++)
    {
     //Weight = TheData(IDataSet,Ipnt,0)*EffN(IDataSet)/TotalN(IDataSet);
     Weight = float(TheData(IDataSet,Ipnt,0));

     // Probability associated observed data (no perfect reader)
     if (Iperfect(IDataSet) == 0)
      {
       Prob1 = 0;
       for (int Age1=MinAge;Age1<=MaxAge;Age1++)
        {
         // Prior probability * product over readers
         Prob2 = Aprob(IDataSet,Age1);
         for (int Ireader=1;Ireader<=Nread(IDataSet);Ireader++)
          if (int(TheData(IDataSet,Ipnt,Ireader)) >= 0)
           {
            AgeA = TheData(IDataSet,Ipnt,Ireader);
            Jreader = ReadPnt(IDataSet,Ireader-1);
            Prob2 *= AgeErrOut(Jreader-1,Age1,AgeA);
           }
         Prob1 += Prob2;
        }
      }

     // Probability associated observed data (a perfect reader)
     if (Iperfect(IDataSet) != 0)
      {
       // True age of this ageing structure
       AgeTrue = TheData(IDataSet,Ipnt,Iperfect(IDataSet));

       // Prior probability * product over readers
       Prob1 = 1;
       for (int Ireader=1;Ireader<=Nread(IDataSet);Ireader++)
        if (TheData(IDataSet,Ipnt,Ireader) >= 0 & Iperfect(IDataSet) != Ireader)
         {
          AgeA = TheData(IDataSet,Ipnt,Ireader);
          Jreader = ReadPnt(IDataSet,Ireader-1);
          Prob1 *= AgeErrOut(Jreader-1,AgeTrue,AgeA);
         }
      }
     Resu(IDataSet,Ipnt) = Prob1;
     Obj_fun(IDataSet) += Weight*log(Prob1+1.0e-100);
    }
   f += Obj_fun(IDataSet);
   //std::cout << Obj_fun(IDataSet) << " " << f << std::endl;
  }

 // Penaly of slopes
 PenalSlope = 0;
 for (int II=0;II<Slope.size();II++) PenalSlope += Slope(II)*Slope(II);

 f = -1*f+AprobWght*Penal + SlopeWght*PenalSlope;

 f += Dummy*Dummy;
 //std::cout << f << std::endl;

 REPORT(TheBias);
 REPORT(TheSD);
 REPORT(Aprob);
 REPORT(AgeErrOut);
 REPORT(f)
 REPORT(Obj_fun);
 return f;
}

