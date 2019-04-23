#include"nanocpp.h"
#include"dft_in_place.h"
#include"tdcmt.h"
using namespace std;


int main(int argc, char* argv[]){

  ////////////////////// INIT THE CODE /////////////////
  nanocpp n(argc,argv);
  cout<<n<<endl; //to screen

  /////////////////// MY CODE ////////////////////
  //init a volume or plane
  //double cds[2*FD_N]={0.5,1.5,0.5,1.5};
  //double cds[2*FD_N]={0.2,1.8,1.5,1.5};
  //int fc[FD_N]={CE,CE};
  //volume<FD_N> kk;
  //kk.init(cds,n.get_grid(),fc);

  /////////////////DFT
  //init a volume
  // double cds[2*FD_N]={0.5,1.5,0.5,1.5};
  // int fc[FD_N]={CE,CE};
  // volume<FD_N> vv;
  // vv.init(cds,n.get_grid(),fc);
  // //init a DFT object
  // int dim=3;
  // double *fw=new double[3];
  // fw[0]=2.75;
  // fw[1]=3.25;
  // fw[2]=4.4;
  // for(int i=0;i<dim;i++)
  //   fw[i]*=2*PI*C0;
  // double dt=n.get_dt();
  // dft_vol df(vv,fw,dim,dt);

  // // Poynting flux with harmonic expansion
  // // Init data
  // double cds[4]={0.,0.99,0.8,0.8};
  // int fc[FD_N]={CE,CE};
  // volume<FD_N> vol;
  // vol.init(cds,n.get_grid(),fc);
  // int nmodes=10;

  //TDCMT
  ifstream fm_in;
  fm_in.open("nmodes.txt", ios::in);
  int tdim[2]; //number of modes they are saved row by row, only even
  char tmp[256];
  for (int i = 0; i != 2; i++){
	fm_in >> tmp;
	tdim[i] = atoi(tmp);
  }
  fm_in.close();  // (x0,z0),(x0,z1),(x0,z2),...
  //reads data from file
  ifstream fp_in;
  fp_in.open("cbd.txt", ios::in);
  double tsizz[4];
  for(int i=0;i<4;i++)
    fp_in>>tsizz[i];
//  double tsizz[6]={0.2,0.802,0.2,0.802, 0.2, 1.804}; //size of vol
  tdcmt td;
  td.setup(tdim,tsizz,n);

  // //HARM TRANSM
  // int dim[1]={1};
  // //double sizz[4]={tsizz[0],tsizz[1],tsizz[3],tsizz[3]};
  // //double sizz[4]={0.,0.99,0.55,0.55};
  // fp_in.open("tfsf.txt", ios::in);
  // for(int i=0;i<4;i++)
  //   fp_in>>tsizz[i];
  // double sizz[4]={tsizz[0],tsizz[1],0.7,0.7};  
  // harm_z ha;
  // ha.init(dim,sizz,n);
  // //HARM REFL
  // //sizz[2]=0.39;//0.35;
  // //sizz[3]=0.39;//0.35;
  // //harm_z ha_r;
  // //ha_r.init(dim,sizz,n);
  fp_in.close();
  

  /////////////////////////////////////////////////////

  ///////////////////// CYCLE /////////////////////
  double sta=MPI_Wtime();
  for(n.tcurr;n.tcurr<=n.total_steps;n.tcurr++){
    
    n.time=n.tcurr*n.dt;     //update time;
    n.single_step();    //single step
    
    /////////////////// MY CODE ////////////////////
    //n.erg_exch(); //exchange boundaries
    //df.accum_dft_all(n,n.tcurr);// accumulate dft results    

    //TDCMT
    td.update_A(n);

    //Harm
    //ha.update_C(n);        
    //ha_r.update_C(n);        

    /////////////////////////////////////////////////////
  }


  /////////////////// MY CODE ////////////////////
  //TDCMT
  string fn=n.get_rootdir()+"/td.bin";
  td.save(fn);
  // //Harm
  // string fn1=n.get_rootdir()+"harm0.bin";
  // ha.save(fn1);
  //fn1=n.get_rootdir()+"harm_r0.bin";
  //ha_r.save(fn1);

  //compute erg and poynting in that region
  //double hh=n.total_erg_in_box(kkk);
  //double hh=n.z_erg_flux(kkk);
  //if(kkk.get_me()==0)
  //  cout<<hh<<endl;
  //  if(kkk.get_me()==0){
  //  fn="res/poy.bin";
  //  poy.save(fn);
  //}
  
  //DFT
  // string fn="res/ergw.bin";
  // df.save_erg(n,fn);
  // fn="res";
  // df.save_fields(n,fn);

  // //fn="res/qy.bin";
  // //n.medium.Qy.save(fn);
  //string fn="res/ez.bin";
  //n.Ez.save<MPI_DOUBLE>(fn,n.get_grid());
  //fn="res/ex.bin";
  //n.Ex.save<MPI_DOUBLE>(fn,n.get_grid());
  //fn="res/hy.bin";
  //n.Hy.save<MPI_DOUBLE>(fn,n.get_grid());
  // //df.save_fields(n,fn);
  
  // delete []fw;
  
  /////////////////////////////////////////////////////  
  
  //////////////////// SAVE OUTPUT AND LOGOUT /////////////

  n.post_process();    //process data
  n.goodbye(sta);   //get time and print goodbye
  
  ////////////////////////////////////////////////////////    


}
