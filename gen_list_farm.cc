void gen_list_farm()
{
  int from = 11;
  int to   = 20;
  int line = 0;
  
  TString str_exe = "";
  TString File = "";
  TString roostr = "";
  
  cout<<endl;
  cout<<endl;

  ofstream ListWrite("cmd_list.csh", ios::out|ios::trunc);
    
  for(int idx=from; idx<=to; idx++)
    {
      line++;
      
      cout<<TString::Format(" ---> processing %5d", line)<<endl;
      
      ofstream WriteFile(TString::Format("./cmd_%06d.csh",idx),ios::out|ios::trunc);
      WriteFile<<"#!/bin/tcsh"<<endl;

      roostr = TString::Format("./read_TFlux -ia %d -ib %d -o %d", (idx-1)*100, idx*100-1, idx);
      
      WriteFile<<str_exe<<roostr<<endl;
      WriteFile<<endl;
      WriteFile.close();

      ListWrite<<TString::Format("nohup ./cmd_%06d.csh > info_%06d &",idx, idx)<<endl;
    }

  ListWrite.close();

  cout<<endl;
  cout<<endl;


}

