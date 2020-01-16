#ifndef __SPACEUTILS__
#define __SPACEUTILS__

///////////////////////////////////////////////////////////////////////

double DeltaFunc(double X, double Y) {
  if(abs(X-Y)<.00001)
    return 1. ;
  else
    return 0. ;
}

///////////////////////////////////////////////////////////////////    

double Gaussian_1D(double mu, double sigma) {
  return exp(-mu * mu /2.0 /sigma /sigma) / sqrt(2.0*M_PI) /sigma ;
}

///////////////////////////////////////////////////////////////////////

double Wrapped_Gaussian(double mu, double sigma, int klim) {
  double sum = 0 ; 
  for(int k=-klim;k<=klim;k++)
    sum += Gaussian_1D(mu-L*(double)k,sigma) ; 
  return sum ;
}

///////////////////////////////////////////////////////////////////////

double Wrapped_Exp(double mu, double sigma, int klim) {
  double sum = 0 ; 
  for(int k=-klim;k<=klim;k++)
    sum += exp( -fabs( L* ( mu+(double)k ) ) / sigma ) ;
  return sum ;
}

///////////////////////////////////////////////////////////////////////

void getCrecCff(char** argv, int nbpop, double *&Crec, double &Cff) {

  Crec = new double [nbpop] ; 
  
  if(IF_CstCrecCff) {
    for(int i=0;i<nbpop;i++) {
      Crec[i] = Sigma[i] ;
      if(IF_Prtr)
	Cff = SigmaExt ;
    }
  }
  else {
    for(int i=0;i<nbpop;i++) 
      Crec[i] = (double) atof(argv[ argIext * nbpop + 6 + i + IF_Prtr ]) ; 
    if(IF_Prtr)
      Cff = (double) atof(argv[ ( argIext + 1) * nbpop + 6 + IF_Prtr ]) ; 
  }
}

///////////////////////////////////////////////////////////////////////

void CreateDir_SpaceCrec(int nbpop,string &path,unsigned long N,double* Crec) {
  
  if(IF_RING) {
    cout  << "O( sqrt(K) ) Cosine interactions " << endl ;
    if(DIM==1)
      if(IF_SPEC)
	path += "/Ring/Spec/" ;
      else
	path += "/Ring/" ;
    else
      path += "/Ring2D/" ;
  }
  if(IF_GAUSS) {
    cout  << "O( sqrt(K) ) Gaussian interactions " << endl ;    
    if(DIM==1)
      path += "/Gauss/" ;
    else
      path += "/Gauss2D/" ;
  }
  if(IF_EXP) {
    cout  << "O( sqrt(K) ) Exp interactions " << endl ;
    if(DIM==1)
      path += "/Exp/" ;
    else
      path += "/Exp2D/" ;
  }

  string mkdirp = "mkdir -p " ;
  
  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;
  string strpop ;
  
  char cCrec[20] ;
  
  if(IF_Dij) {
    for(int i=0;i<nbpop;i++)    
      for(int j=0;j<nbpop;j++) {
	sprintf(cCrec,"%0.4f",Crec[j] * Dij[j+i*nbpop]) ; 
	printf(cCrec,"%0.4f ",Crec[j] * Dij[j+i*nbpop]) ; 	
	path += "Crec"+ popList[i] + popList[j] + string(cCrec) ; 
      }
  }
  else {
    for(int i=0;i<nbpop;i++) {
      strpop = popList[i] ;
      sprintf(cCrec,"%0.4f",Crec[i]) ; 
      path += "Crec"+ strpop + string(cCrec) ; 
    }    
  }
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str() ;
  const int dir_err = system(cmd) ; 

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////////

void CreateDir_SpaceCrec2D(int nbpop,string &path,int N,double* Crec) {
  
  if(IF_RING) {
    cout  << "O( sqrt(K) ) Cosine interactions " << endl ;
    path += "/Ring2D/" ;
  }
  if(IF_GAUSS) {
    cout  << "O( sqrt(K) ) Gaussian interactions " << endl ;
    path += "/Gauss2D/" ;
  }
  string mkdirp = "mkdir -p " ;
  
  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;
  string strpop ;
  
  char cCrec[10] ;
  
  for(int i=0;i<nbpop;i++) {
    strpop = popList[i] ;
    sprintf(cCrec,"%0.2f",Crec[i]) ; 
    path += "Crec"+ strpop + string(cCrec) ; 
  }
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////////

void CreateDir_SpaceCrecCab(int nbpop,string &path,int N,double **Crec) {
  
  if(IF_RING) {
    cout  << "O( sqrt(K) ) Cosine interactions " << endl ;
    path += "/Ring/" ;
  }
  if(IF_GAUSS) {
    cout  << "O( sqrt(K) ) Gaussian interactions " << endl ;
    path += "/Gauss/" ;
  }
  if(IF_EXP) {
    cout  << "O( sqrt(K) ) Exp interactions " << endl ;
    path += "/Exp/" ;
  }

  string mkdirp = "mkdir -p " ; 
  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;
  string strPre,strPost ;

  char cCrec[10] ;

  for(int i=0;i<nbpop;i++) {
    strPost = popList[i] ;
    // cout << strPost << " " << strPre  << endl ;
    for(int j=0;j<nbpop;j++) {
      strPre = popList[j] ;
      sprintf(cCrec,"%0.4f",Crec[i][j]) ;
      path += "Crec"+ strPost + strPre + string(cCrec) ; 
    }
  }
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

void CreateDir_SpaceCff(int nbpop,string &path,int N,double Cff) {
  
  string mkdirp = "mkdir -p " ;
 
  char cCff[10] ;
  sprintf(cCff,"%0.4f",Cff) ;
  string sCff = string(cCff) ;
  
  path += "Cff"+sCff ; 
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////////

void External_Input(int nbpop, unsigned long N, unsigned long* Cpt, double K, double Cff, double* Iext, double *IextBL, double *&IextFF, string path) { 
/* void External_Input(int nbpop, unsigned long N, unsigned long* nbN, double K, double Cff, double* Iext, double *IextBL, int PrtrPop, vector<vector<double> > &Jext, string path) { */

  random_device rd ; 
  default_random_engine gen( rd() ) ; 
  uniform_real_distribution<double> unifExt(0, 1.0) ; 
  normal_distribution<double> gaussianExt(0, SIGMA_EXT ) ; 
  
  cout << "External Input : " ; 

  double *X ; 
  X = (double *) malloc( (unsigned long) N * sizeof(double)) ;
  
  for(int i=0;i<nbpop;i++) 
    for(unsigned long j=Cpt[i];j<Cpt[i+1];j++) 
      X[j] = L * fmod( double(j-Cpt[i]), double(Cpt[i+1]-Cpt[i]) ) / ( double(Cpt[i+1]-Cpt[i]) ) ; 

  for(int i=0;i<nbpop;i++) 
    cout << X[Cpt[i]] << " " << X[(Cpt[i+1]-Cpt[i])/2] << " " << X[Cpt[i+1]-1] << " "; 
  cout << endl ; 

  double *p ; 
  p = (double *) malloc( (double) nbpop * sizeof(double)) ;

  int dum = 0 ;
  for(int i=0;i<nbpop;i++) {
    p[i] = Iext[i] - IextBL[i] ; 
    if( p[i]>.000000001 || p[i]<0 ) 
      dum = 1 ; 
    if(IF_WEAK) 
      p[i] = p[i]/sqrt(K) ; 
  } 
  
  cout << "Baseline " ;
  for(int i=0;i<nbpop;i++) 
    cout << IextBL[i]/sqrt(K)/m0/(Vth-Vr) << " " ; 
  
  if(dum==0) { 
    if(PrtrProfile==-1) { 
      /* cout << "| FeedForward On " << endl ;  */
      /* for(int i=0;i<nbpop;i++)  */
      /* 	for(unsigned long j=Cpt[i];j<Cpt[i+1];j++)  */
      /* 	  IextFF[j] = IextBL[i] * Cff * cos( 2. * ( X[j] - X[(Cpt[i+1]-Cpt[i])/2-1] - PHI0 ) ) ;  */
    } 
    else 
      cout << endl ; 
  }
  else { 
    
    cout << "| Perturbation " ; 
    for(int i=0;i<nbpop;i++) 
      cout << Iext[i]/sqrt(K)/m0/(Vth-Vr) - IextBL[i]/sqrt(K)/m0/(Vth-Vr) << " " ; 
    cout << endl ; 

    if(PrtrPop>=0) {
      for(unsigned long j=Cpt[PrtrPop];j<Cpt[PrtrPop+1];j++) {
	switch(PrtrProfile) { 
	case -1 :
  	  IextFF[j] = Iext[PrtrPop] * Cff * cos( 2. * ( X[j] - L/2 ) ) ; 

	case 0 : 
	  if( abs( X[j]-L/2 ) <= Cff * L / 2.0 ) 
	    IextFF[j] = p[PrtrPop]*(1.0-IF_OPSIN) + p[PrtrPop]*IF_OPSIN*gaussianExt(gen) ; 	
	  break ; 
	  
	case 1 :
	  // IextFF[j] = p[PrtrPop]*Gaussian_1D(X[j]-L/2,Cff) ; 
	  IextFF[j] = p[PrtrPop]*Wrapped_Gaussian(X[j]-L/2,Cff,4) ; 
	  if(IF_LIGHT_SCATTER)
	    IextFF[j] = p[PrtrPop]*( Wrapped_Gaussian(X[j]-L/2,Cff,4) + 
				     Wrapped_Gaussian(X[j]-L/2,SIG_SCATTER,4) / sqrt(K) ) ; 

	  break ; 
	
	case 2 : 
	  IextFF[j] = p[PrtrPop]*exp(-abs(X[j]-L/2)/Cff)/Cff ; 
	  break ; 

	default : 
	  if( abs(X[j]-X[(Cpt[PrtrPop+1]-Cpt[PrtrPop])/2-1]) <= Cff / 2.0 ) 
	    IextFF[j] = p[PrtrPop] ; 
	  break ; 
	} 
      } 

      /* for(unsigned long j=Cpt[PrtrPop+2];j<Cpt[PrtrPop+1+2];j++) { */
      /* 	switch(PrtrProfile) { */
      /* 	case 0 : */
      /* 	  if( abs( X[j]-L/2 ) <= Cff * L / 2.0 ) */
      /* 	    IextFF[j] = p[PrtrPop]*(1.0-IF_OPSIN) + p[PrtrPop]*IF_OPSIN*gaussianExt(gen) ; */
      /* 	  break ; */
	  
      /* 	case 1 : */
      /* 	  // IextFF[j] = p[PrtrPop]*Gaussian_1D(X[j]-L/2,Cff) ; */
      /* 	  IextFF[j] = p[PrtrPop]*Wrapped_Gaussian(X[j]-L/2,Cff*2.0,4) ; */
      /* 	  break ; */
	  
      /* 	default : */
      /* 	  if( abs(X[j]-X[(Cpt[PrtrPop+1]-Cpt[PrtrPop])/2-1]) <= Cff / 2.0 ) */
      /* 	    IextFF[j] = p[PrtrPop] ; */
      /* 	  break ; */
      /* 	} */
      /* } */
      
    }
    else {
      for(int i=0;i<nbpop;i++) 
	for(unsigned long j=Cpt[i];j<Cpt[i+1];j++) {
	  switch(PrtrProfile) { 
	  case 0 : 
	    if( abs( X[j]-L/2 ) <= Cff * L / 2.0 ) 
	      IextFF[j] = p[i]*(1.0-IF_OPSIN) + p[i]*IF_OPSIN*gaussianExt(gen) ; 	
	    break ; 
	    
	  case 1 :
	    // IextFF[j] = p[i]*Gaussian_1D(X[j]-L/2,Cff) ; 
	    IextFF[j] = p[i]*Wrapped_Gaussian(X[j]-L/2,Cff,4) ; 
	    break ;
	    
	  default :
	    if( abs(X[j]-X[(Cpt[i+1]-Cpt[i])/2-1]) <= Cff / 2.0 ) 
	      IextFF[j] = p[i] ;
	    break ;
	  }
	}
    }
  }
  
  free(X) ; 
  free(p) ;
}

///////////////////////////////////////////////////////////////////////

void External_Input2D(int nbpop, unsigned long N, unsigned long* nbN, double K, double Cff, double* Iext, double *IextBL, int ndI, double *&Jext, string path) {

  cout << "External Input : " ;

  /* double **X ; */
  /* X = new double*[nbpop] ;       */
  /* for(int i=0;i<nbpop;i++) */
  /*   X[i] = new double[nbN[i]] ; */

  /* for(int i=0;i<nbpop;i++) */
  /*   for(unsigned long j=0;j<nbN[i];j++)  */
  /*     X[i][j] = L * fmod( double(j), sqrt( double( nbN[i]) ) ) / sqrt( double( nbN[i] ) ) ; */
  
  /* double **Y ; */
  /* Y = new double*[nbpop] ;       */
  /* for(int i=0;i<nbpop;i++) */
  /*   Y[i] = new double[nbN[i]] ; */

  /* for(int i=0;i<nbpop;i++) */
  /*   for(unsigned long j=0;j<nbN[i];j++) */
  /*     Y[i][j] = L * floor( double(j) / sqrt( double( nbN[i]) ) ) / sqrt( double (nbN[i]) ) ; */
  
  /* double p = 1.0 ; */
  
  /* cout << "Baseline " ; */
  /* for(int i=0;i<nbpop;i++) */
  /*   cout << IextBL[i]/sqrt(K)/m0 << " " ; */
  
  /* if(abs(Iext[ndI]-IextBL[ndI])<=.001) { */
  /*   cout << endl ; */
  /* } */
  /* else { */
  /*   p = Iext[ndI] - IextBL[ndI] ; */
    
  /*   cout << "| Perturbation " ; */
  /*   for(int i=0;i<nbpop;i++) */
  /*     cout << Iext[i]/sqrt(K)/m0 - IextBL[i]/sqrt(K)/m0 << " " ; */
  /*   cout << endl ; */

  /*   for(unsigned long j=0;j<nbN[ndI];j++) { */

  /*     switch(IF_GAUSS) { */
	
      /* case 0 : */
      /* 	if( ( X[ndI][j]-.5*L ) * ( X[ndI][j]-.5*L ) + ( Y[ndI][j]-.5*L ) * ( Y[ndI][j]-.5*L ) <= Cff*Cff / 4.0 )  */
      /* 	  Jext[ndI][j] = p ;  */
      /* 	else */
      /* 	  Jext[ndI][j] = 0. ;  */
      /* 	break ; */
      
      /* /\* case 1 : *\/ */
	
      /* /\* 	if( abs( X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= 4.0 * Cff ) *\/ */
      /* /\* 	    Jext[ndI][j] = p*PeriodicGaussian(X[ndI][j],X[ndI][nbN[ndI]/2-1],Cff) ; *\/ */
      /* /\* 	  else *\/ */
      /* /\* 	    Jext[ndI][j] = 0. ; *\/ */
      /* /\* 	break ; *\/ */

      /* /\* case 2 : *\/ */

      /* /\* 	if( abs(X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= Cff ) *\/ */
      /* /\* 	  Jext[ndI][j] = p ; // half height of the Gaussian *\/ */
      /* /\* 	else *\/ */
      /* /\* 	  if( abs(X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= 4.0 * Cff ) *\/ */
      /* /\* 	    Jext[ndI][j] = p*PeriodicGaussian(X[ndI][j],X[ndI][nbN[ndI]/2-1],Cff) / ( exp(- 1.0 / 2.0) / sqrt(2.0*M_PI) / Cff ) ; *\/ */
      /* /\* 	  else *\/ */
      /* /\* 	    Jext[ndI][j] = 0. ; *\/ */
      /* /\* 	break ; *\/ */

      /* default : */
      /* 	if( abs( X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= Cff / 2.0 )  */
      /* 	  Jext[ndI][j] = p ;  */
      /* 	else */
      /* 	  Jext[ndI][j] = 0. ;  */

      /* 	break ; */
	
  /*     }	 */
      
  /*   } */
    
  /* } */
        
  /* delete [] X ; */
  /* delete [] Y ; */
}

//////////////////////////////////////////////////////////////////

#endif
