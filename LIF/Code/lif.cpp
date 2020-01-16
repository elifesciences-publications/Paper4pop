#include "librairies.h"
#include "GlobalVars.h" 
#include "Net_Utils.h"
#include "Matrix_Utils.h"
#include "Space_Utils.h"
#include "MFrates.h"

clock_t t1=clock();

int main(int argc , char** argv) {

  string dir ;
  int nbpop, ndI ; 
  unsigned long N ;  
  double dt, K, **J, **Tsyn, g, *Iext, *IextBL, *Crec, Cff = 4.0*L ; 
  double *RatesMF ;

  bool IF_STRUCTURE=0 ;  
  IF_STRUCTURE = IF_RING || IF_GAUSS || IF_EXP ; 

  dt = DT ; 
  if(IF_BENCHMARK) 
    dt = (double) atof(argv[6]) ;

  Set_Parameters(argc,argv,dir,nbpop,N,K,g,Iext,IextBL) ; 

  if(IF_Prtr) { 
    ndI = min(nbpop-1,PrtrPop) ; 
    if(ndI>=0)
      Iext[ndI] += (double) atof(argv[6]) ; 
    else {
      for(int i=0;i<nbpop;i++) 
	Iext[i] += (double) atof(argv[6]) ; 
      ndI = 0 ;
    }

    if(IF_PROP) 
      Cff = (double) atof(argv[7]) ;     
  }

  cout << "Synaptic Strength : " << endl ; 
  Import_Connectivity_Parameters(nbpop,J,dir) ; 
  for(int i=0;i<nbpop;i++) { 
    for(int j=0;j<nbpop;j++) 
      cout << J[i][j] << " "; 
    cout << endl ; 
  } 

  MF_Rates(nbpop,IextBL,J,RatesMF) ;

  cout << "MF Rates: " ; 
  for(int i=0;i<nbpop;i++) 
    cout << RatesMF[i] << " " ; 
  cout << endl ; 

  ///////////////////////////////////////////////////////////////////    
  // Time Csts 
  ///////////////////////////////////////////////////////////////////    
  
  cout << "Membrane time constants : " ;
  for(int i=0;i<nbpop;i++) 
    cout << Tm[i] << " " ;  
  cout << endl ;    
    
  Import_Synaptic_Parameters(nbpop,Tsyn,dir) ; 

  double **CstTsyn = new double *[nbpop]() ;
  for(int i=0;i<nbpop;i++) {
    CstTsyn[i] = new double [nbpop]() ; 
    for(int j=0;j<nbpop;j++) 
      CstTsyn[i][j] = exp(-dt/Tsyn[i][j]) ; 
  }

  int i=0,j=0 ; 
  unsigned long k=0, l=0 ; 
  
  ///////////////////////////////////////////////////////////////////  
  // Connectivity Matrix
  ///////////////////////////////////////////////////////////////////  

  string Jpath = "../../" ;
  Create_Path(nbpop,Jpath,N,K) ;   

  if(IF_STRUCTURE) {
    getCrecCff(argv,nbpop,Crec,Cff) ; 
    CreateDir_SpaceCrec(nbpop,Jpath,N,Crec) ; 
  }
  
  
  int *nbPost ;
  unsigned long *idxPost ;
  unsigned long *IdPost ; 
  Import_Connectivity_Matrix(nbpop,N,Jpath,nbPost,idxPost,IdPost) ; 

  ///////////////////////////////////////////////////////////////////
  // Path 
  ///////////////////////////////////////////////////////////////////    

  string path = "../" ; 
  CreateDir(dir, nbpop, N, K, g, path) ; 

  if(IF_BENCHMARK) {

    char cdt[10] ;
    sprintf(cdt,"%0.4f",dt) ; 
    string sdt = string(cdt) ; 
    
    path += "/BenchMark/dt" + sdt ; 

    string mkdirp = "mkdir -p " ; 
    mkdirp += path ; 

    const char * cmd = mkdirp.c_str() ; 
    const int dir_err = system(cmd) ; 
    
    if(-1 == dir_err) 
      cout << "error creating directories" << endl ; 
    cout << path << endl ;    
  }

  if(IF_SHARED)
    CreateDirShared(path) ;

  if(IF_STRUCTURE) { 
    CreateDir_SpaceCrec(nbpop,path,N,Crec) ; 
    if(IF_Prtr) 
      CreateDir_SpaceCff(nbpop,path,N,Cff) ; 
     
  }
  string popList[4] = {"E","I","S","V"} ; 
  if(nbpop==1) 
    popList[0] = "I" ; 

  if(IF_Prtr) { 
    if(PrtrPop>=0)
      cout <<"Perturbed Pop " << popList[ndI] << " | Input "<< Iext[ndI]-IextBL[ndI] << endl ; 
    if(IF_PROP) {
      cout << Cff << " only "  << endl ; 
      CreateDir_Prop(nbpop,ndI,Iext[ndI],Cff,path) ; 
    } 
    else 
      CreateDir_Iext(nbpop,ndI,Iext[ndI],path) ; 
  } 
  
  if(IF_JabLoop) { 

    if(IF_M0) { 
      
    } 
    else {
       Iext[0] = Je0 ; 
       IextBL[0] = Je0 ; 
    }
    
    J[Ax][Ay] = (double) atof(argv[ IF_STRUCTURE * (nbpop + 1) + 6 + IF_Prtr ]) ; 
    cout <<" J_ " << popList[Ax] << popList[Ay] << " " << J[Ax][Ay] << endl ;
    CreateDir_JabLoop(path, J[Ax][Ay]) ; 
  }

  ///////////////////////////////////////////////////////////////////    
  // number of neurons
  ///////////////////////////////////////////////////////////////////    

  unsigned long *nbN ;
  nbNeurons(nbpop,N,nbN) ;
  unsigned long *Cpt ; // compteur Cpt0 = Ne, cpt1 = Ne+Ni ...
  cptNeurons(nbpop,nbN,Cpt) ;

  Save_Parameters(dir, nbpop, N, nbN, dt, K, J, Iext, Tm, Tsyn, path) ;
   
  int *whichPop ;
  whichPop = (int *) malloc( (unsigned long long) N * sizeof(int) ) ;

  for(i=0;i<nbpop;i++) 
    for(k=0;k<N;k++)
      if(k>=Cpt[i] && k<Cpt[i+1])
	whichPop[k] = i ; 
  
  ///////////////////////////////////////////////////////////////////     
  // Scaling
  //////////////////////////////////////////////////////////////////

  for(i=0;i<nbpop;i++) {
    Iext[i] = sqrt(K) * Iext[i] * m0 * (Vth-Vr) ; 
    IextBL[i] = sqrt(K) * IextBL[i] * m0 * (Vth-Vr) ; 
    for(j=0;j<nbpop;j++) 
      J[i][j] = g * J[i][j] / sqrt(K) / Tsyn[i][j] * (Vth-Vr) ; 
  } 

  if(IF_EtoSOM_O1 && nbpop>2) J[2][0] = J[2][0] / sqrt(K) ; 
  if(IF_PROP) 
    if(IF_WEAKPROP)
      Cff = Cff/sqrt(K) ; 

  ///////////////////////////////////////////////////////////////////    
  // Allocating variables
  ///////////////////////////////////////////////////////////////////    

  vector<vector<vector<double> > >Isyn(nbpop,vector<vector<double> >(nbpop)) ; 

  for(i=0;i<nbpop;i++)
    for(j=0;j<nbpop;j++) 
      Isyn[i][j].resize(nbN[i]) ; 
  
  double *Mean_Activity ; 
  double *Idv_Activity, *Itot, *ItotRK2 ;
  double *IsynDecay, *IsynDecayUpd, *IsynRise, *IsynRiseUpd ; 
  double *Volt, *Vold, *tspk ;
  double *IextFF ; 

  Mean_Activity = (double *) malloc( (int) nbpop * sizeof(double)); 
  Idv_Activity = (double *) malloc( (unsigned long long) N * sizeof(double)); 
  Itot = (double *) malloc( (unsigned long long) N * sizeof(double)); 
  ItotRK2 = (double *) malloc( (unsigned long long) N * sizeof(double)); 
  IextFF = (double *) malloc( (unsigned long long) N * sizeof(double)); 

  Volt = (double *) malloc( (unsigned long long) N * sizeof(double)); 
  Vold = (double *) malloc( (unsigned long long) N * sizeof(double)); 
  tspk = (double *) malloc( (unsigned long long) N * sizeof(double)); 

  ///////////////////////////////////////////////////////////////////    
  // opening files 
  ///////////////////////////////////////////////////////////////////    

  string strMean = path + "/Mean.txt" ; 
  ofstream fMean(strMean.c_str(), ios::out | ios::ate);
 
  string strIdvRates = path + "/IdvRates.txt" ; 
  ofstream fIdvRates(strIdvRates.c_str(), ios::out | ios::ate);
 
  string strVoltage = path + "/Voltage.txt" ; 
  ofstream fVoltage(strVoltage.c_str(), ios::out | ios::ate);

  string strRaster = path + "/Raster.txt" ; 
  ofstream fRaster(strRaster.c_str(), ios::out | ios::ate);
  
  ///////////////////////////////////////////////////////////////////     
  cout << "Initialization" << endl;
  ///////////////////////////////////////////////////////////////////    
 
  random_device rd ; 
  default_random_engine gen( rd() ) ; 
  uniform_real_distribution<double> unif(0,1) ; 
  normal_distribution<double> gaussianI( 0, (Vth-Vr)/4.0 ) ; 
  uniform_real_distribution<double> unifV(Vr,Vth) ; 

  cout << "Check random seed " ; 
  for(i=0;i<10;i++) 
    cout << unif(gen) << " " ; 
  cout << endl ; 
  
  for(k=0;k<N;k++) {
    Volt[k] = Tm[whichPop[k]]*Vr ; 
    Itot[k] = IextBL[whichPop[k]] ; 
    ItotRK2[k] = Itot[k] ; 
    IextFF[k] = 0.0 ; 
  }

  for(i=0;i<nbpop;i++) 
    for(j=0;j<nbpop;j++)
      for(k=0;k<nbN[i];k++) { 
	Isyn[i][j][k] = (K*J[i][j]*Tsyn[i][j]*RatesMF[j]  / (double) nbN[i] + gaussianI(gen) /Tm[i] ) ; 
	Itot[k+Cpt[i]] += Isyn[i][j][k] ; 
	ItotRK2[k+Cpt[i]] += CstTsyn[i][j]*Isyn[i][j][k] ; 
      }
  
  ///////////////////////////////////////////////////////////////////     
  // External Perturbation
  ///////////////////////////////////////////////////////////////////     

  if(IF_Prtr) 
    if(DIM==1) 
      External_Input(nbpop,N,Cpt,K,Cff,Iext,IextBL,IextFF,path) ; 
    else 
      External_Input2D(nbpop,N,nbN,K,Cff,Iext,IextBL,ndI,IextFF,path) ; 
  
  /////////////////////////////////////////////////////////////////// 

  cout << "Scheme: " ; 
  if(IF_EULER) cout << "EULER " ; 
  if(IF_RK2) cout << "RK2 " ; 
  if(IF_INTERPOLATION) cout << "with interpolation " ; 
  cout << "and exponential decay" ;
  cout << endl ; 
  
  /////////////////////////////////////////////////////////////////// 
  // Main Loop 
  ///////////////////////////////////////////////////////////////////     

  double tw=0., tc=0. ; //Time window 
  double RK1=0,RK2=0 ;
  double percentage=0 ; 
  double tTRANS=0 ;
  int IF_TRANSIENT = 0 ;

  if(IF_TRANSIENT_IEXT) {
    IextBL[2] = -100.0*sqrt(K)*m0*(Vth-Vr) ; 
    cout << "Transient Input" << endl ;
  }

  int PreS=0, Post=0 ; 
  
  cout << "Main loop :" ; 
  cout << " duration " << duration << " | dt " << dt ;
  cout << " | Tst " << Tst << " | Tw " << Tw << " | Tl " << Tl << endl ; 

  for (double t=0.; t<=duration+Tst+dt; t+=dt) { 

    percentage = t/(duration+Tst) ;

    if(IF_TRANSIENT_IEXT && t>Tst/2.0) 
      while(IextBL[2]<=0.0) 
	IextBL[2] += dt*100*sqrt(K)*m0*(Vth-Vr)/nbIext ; 

    if(t>=duration+Tst-Tl ) fVoltage << t-duration-Tst+Tl ; // Adding Spikes by hand 
    
    //Updating Voltage 
    for (k=0; k<N; k++) { 
      PreS = whichPop[k] ;
      Vold[k] = Volt[k] ; // V(t) 

      if(IF_EULER) // Standard Euler Algorithm 
	Volt[k] = (1.0-dt/Tm[PreS]) * Volt[k] + dt * ( Itot[k] + Vr / Tm[PreS] ) ; // V(t+dt) 
      
      if(IF_RK2) { // RK2 	  	
	RK1 = -(Volt[k]-Vr)/Tm[PreS] + Itot[k] ; 
	RK2 = -(Volt[k]-Vr+dt*RK1)/Tm[PreS] + ItotRK2[k] ; 
	Volt[k] = Volt[k] + dt/2.0 * ( RK1 + RK2 ) ; 
      } 
      
      // reset total input
      Itot[k] = IextBL[PreS] + IextFF[k] ; 
      ItotRK2[k] = Itot[k] ; 

      if (Volt[k]>=Vth) { //if Spike reset V and update postsynaptic currents 
	
	if(IF_INTERPOLATION) { // Improved accuracy 
	  tspk[k] = t + dt * (Vth-Vold[k]) / (Volt[k]-Vold[k]) ; 
	  Volt[k] = (Volt[k]-Vth) * (1.0 + dt/Tm[PreS]*(Vold[k]-Vr) /(Volt[k]-Vold[k]) ) + Vr ; 
	} 
	else { 
	  tspk[k] = t ; 
	  Volt[k] = Vr ; 
	} 
	
	if(t>=Tst) { 
	  Mean_Activity[PreS] += 1. ; 
	  Idv_Activity[k] += 1. ; 
	  if(k<10 && t>=duration+Tst-Tl) fVoltage << " " << Vpeak ; 
	} 

	// Updating Postsynaptics inputs 
	if(IF_INTERPOLATION) { 
	  for(l=idxPost[k];l<(unsigned long)idxPost[k]+(unsigned long)nbPost[k];l++) { 
	    Post = whichPop[IdPost[l]] ; 
	    if(J[Post][PreS]!=0) 
	      Isyn[Post][PreS][IdPost[l]-Cpt[Post]] += J[Post][PreS]*exp(-(t-tspk[k])/Tsyn[Post][PreS]) ;
	  } 
	} 
	else { 
	  for(l=idxPost[k];l<(unsigned long)idxPost[k]+(unsigned long)nbPost[k];l++) { 
	    Post = whichPop[IdPost[l]] ; 
	    if(J[Post][PreS]!=0) 
	      Isyn[Post][PreS][IdPost[l]-Cpt[Post]] += J[Post][PreS] ; 
	  } 
	} 
	// Updating ISI 
	if(t>=duration+Tst-Tl) 
	  fRaster << fixed << setprecision(1) << (float) (k) << " " << (float) (tspk[k]-duration-Tst+Tl) << endl ; 
	
      } // endif spike 
      else 
	if(k<10 && t>=duration+Tst-Tl) fVoltage << " " << Volt[k] ; 
    } // endfor neurons 
    
    // Updating Total Synaptic Input to each neurons in each population 
    for(i=0;i<nbpop;i++)
      for(j=0;j<nbpop;j++) 
	for(k=0;k<nbN[i];k++) {
	  Isyn[i][j][k] *= CstTsyn[i][j] ; 
	  Itot[k+Cpt[i]] += Isyn[i][j][k] ; 
	  ItotRK2[k+Cpt[i]] += CstTsyn[i][j]*Isyn[i][j][k] ; // Isyn(t+dt) as if no spike emitted 
	}
    
    if(t>=duration+Tst-Tl) fVoltage << endl ; 
    
    if(tw>=Tw) {
      cout << int(percentage*100.0) << "% " ; 
      cout << "t " << int(t-Tst-Tw) << " Rates" ; 
      for(i=0;i<nbpop;i++) 
    	cout << " " << Mean_Activity[i] / tw*1000. / (double) nbN[i] ; 
      cout << "\r"  ;
      cout.flush() ; 
        
      fMean << t-Tst-Tw ; 
      for(i=0;i<nbpop;i++) {
    	fMean << " " << Mean_Activity[i] / tw*1000. / (double) nbN[i] ; 
	Mean_Activity[i] = 0 ; 
      }
      fMean << endl ; 
    
      fIdvRates << t-Tst-Tw ; 
      for(k=0 ;k<N;k++) { 
    	fIdvRates <<" "<< Idv_Activity[k] / tw*1000. ; 
    	Idv_Activity[k] = 0 ; 
      } 
      fIdvRates << endl ; 
      
      tw=0 ;
      
    }//ENDIF 
    else {
      cout << int(percentage*100.0) << "% " << "\r" ;
      cout.flush() ;
    }
    // updating Time Windows 
    if(t>=Tst) {
      tw += dt ;
      if(IF_TIMECOURSE) {
	tc += dt ;
	
	if(tc>=Tc) {
	  if(tc>Tc+dTc+Tw) { 
	    if(IF_TRANSIENT) {
	      cout << endl << "Feedforward Off " << endl ; 	      
	      IextBL[PrtrPop] -= sqrt(K) * dPrtr * m0 * (Vth-Vr) ;
	      IF_TRANSIENT = 0 ; 
	    }
	  }
	  else {
	    if(!IF_TRANSIENT) {	
	      IextBL[PrtrPop] += sqrt(K) * dPrtr * m0 * (Vth-Vr) ;
	      cout << endl << "Feedforward On " << PrtrPop << " " << IextBL[PrtrPop] << endl ;
	      IF_TRANSIENT = 1 ;
	    }  
	  } 
	}
      }

    } 
  } //ENDMAINLOOP 

  cout << endl ; 

  ///////////////////////////////////////////////////////////////////    

  Isyn.clear() ;

  free(Itot) ; 
  free(ItotRK2) ;  
  free(IextFF) ;

  free(Volt) ; 
  free(Vold) ; 
  free(tspk) ;

  free(Mean_Activity) ;
  free(Idv_Activity) ; 

  ///////////////////////////////////////////////////////////////////    

  fMean.close(); 
  fIdvRates.close();
  fVoltage.close();
  fRaster.close();

  ///////////////////////////////////////////////////////////////////    

  cout << "Simulation Done !" << endl ; 

  clock_t t2=clock();  
  int HOURS=0,MIN=0,SEC=0;
  string str_TIME = path + "/CPU_TIME.txt" ; 
  ofstream TIME(str_TIME.c_str(), ios::out | ios::ate);

  SEC = (t2-t1)/CLOCKS_PER_SEC ;
  HOURS = SEC/3600 ;
  MIN = SEC/60 ;
  SEC = SEC % 60 ;
  cout << "Elapsed Time = " << HOURS << "h " << MIN << "m " << SEC << "s" << endl;
  TIME << "Elapsed Time = " << HOURS << "h " << MIN << "m " << SEC << "s" << endl;
  TIME.close() ;

  return 0;

}
