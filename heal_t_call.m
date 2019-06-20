function heal_t_call(inppath,outpath,sel,selp,fss,ini_filt,inc_filt,overlap,W_num,s_sel)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% HEAL-T is designed for users who are comfortable using the MATLAB 
% environment to run software but does not require advanced programing 
% knowledge. 
%
% Authors: 
% Juan Manuel Mayor Torres (juan.mayortorres@unitn.it)
%
% If the user chooses to use the HEAL-T to clean and debug PPG signal
% channels in clinical and in-vivo trials the user should cite the following paper:
%
% Torres, J. M. M., Ghosh, A., Stepanov, E. A., and Riccardi, G. (2016). HEAL-T: An efficient PPG-based heart-rate and IBI estimation method during physical exercise. 
% In Signal Processing Conference (EUSIPCO), 2016 24th European (pp. 1438-1442). IEEE
% 
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% inppath: is the input path where the bvp files,folder_pathccel and/or GSR data are allocated.
%% outpath: outpath given by default or by the user
%% sel: ibifolder selection based on the method 0-> IBIHEALPEAK 1-> ibi.txt
%% sel_p: selector for HEAL-T method and the baseline, 1 -> HEAl-T method, 0 -> baseline
%% fss: please specify the sampling frequency of your data (Hz).
%% ini_filt: initial BHW bandwidth (default [0.7 2.5]);
%% inc_filt: BHW constant increment value (default 0);
%% overlap: For old HEAL-T version: percentage of overlap (e.g 0.5 , default 1 or no overlap). For the new version of HEAL-T: this number is the overlap in seconds, in general this is the length of the overlap window.
%% W_num: For old HEAL-T version: the number of windows for segment the BVP signal (e.g 50). For the new version of HEAL-T: this number is effective length of the window in seconds.
%% s_sel: spline type selector 1-> Cubic Spline, 0-> Moving Average Filter
%% Initial Recursive bvp.txt and ibi.txt file searching..

addpath(genpath([pwd '/functions_HR/thirdparty']));

if (nargin==5)
   ini_filt{1}=[0.7 2.5];
   inc_filt=[0 0];
elseif (nargin<5)
    disp('Not enough inputs!!');   
end;

display(inppath)
av_in=1;
av_k=0;
folder_path=dir(inppath);
if (length(folder_path)==0)
    disp(['Folder Path does not exist!!..Please specify another location']);
end;
for (i=1:length(folder_path))
  if (folder_path(i).isdir==1 && all(folder_path(i).name~='.')) 
    heal_t_call([inppath '/' folder_path(i).name],outpath,sel,selp,ini_filt);
  else
    for (j=1:length(folder_path))
        if (selp==0)
            if(sel==0)
                if (length(folder_path(j).name)==7 && all(folder_path(j).name=='ibi.txt'))
                    if (av_k==0)
                        av_in=0;
                        disp('File Exist!');
                    else
                        av_in=1;
                    end;
                    break;
                end;
            else
                if (length(folder_path(j).name)==11 && all(folder_path(j).name=='iibibas.txt'))
                    av_in=0;
                    disp('File Exist!');
                    break;
                end;
               end;
        else
            if (length(folder_path(j).name)==16 && all(folder_path(j).name=='ibiIBIselHEB.txt'))
                av_in=0;
                disp('File Exist!');
                break;
            end;  
        end;
   end;
%% asking the user for a results update
   if (av_in==0)
    ind_val=input('Do you want to update your current output file? [Y/N]:','s');
    av_in=1;
     if (ind_val~='y' || ind_val~='Y')
      temp='bvp.txt';
     else
        break;
     end;
   else
       temp=folder_path(i).name;
   end;
   if (length(temp)==7 && all(temp=='bvp.txt') && av_in==1)
          pos_slash=find(inppath=='/');
          if (selp==0)
           k=0;
           read_file_HEAL_no_over(inppath(1:pos_slash(length(pos_slash)-3)-1),outpath,inppath(pos_slash(length(pos_slash)-3)+1:pos_slash(length(pos_slash)-1)-1),inppath(pos_slash(length(pos_slash)-1)+1:length(inppath)),{'bvp','acc'},10,W_num,0,sel,fss,ini_filt,inc_filt,overlap,s_sel);
           break;
           av_k=1;
           %% this is for run the old version of IBI and HR detection based on single filter criterion and not a robust HR spectral peak searching (old version***)
          else
           k=1;
           runsubjectnew_heal_t(inppath(1:pos_slash(length(pos_slash)-3)-1),outpath,inppath(pos_slash(length(pos_slash)-3)+1:pos_slash(length(pos_slash)-1)-1),inppath(pos_slash(length(pos_slash)-1)+1:length(inppath)),{'bvp','acc'},fss,ini_filt,inc_filt,W_num,overlap);
           break;
           %% this is new method based on a robust HR searching process and dynamic BHW filter's bandwidth change (new HEAL-T*).
        end; 
     else
      %disp('bvp.txt File does not exist!!')
    end;
   end
end;

