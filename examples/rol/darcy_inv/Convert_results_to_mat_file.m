clear
close all
clc

vary_time_steps = true;
vary_rank_and_storage = false;

if vary_time_steps

    mkdir vary_time_steps

    params = Read_Inputs_from_Xml();
    if params.N_fe > 0
        fe_results = Read_Results('Forward_Euler_Cost_Report.txt');
        filename = ['FE_with_',num2str(params.N_fe),'_time_steps.mat'];
        cd vary_time_steps
        save(filename,'params','fe_results');
        cd ../
    end

    if params.N_me > 0
        me_results = Read_Results('Modified_Euler_Cost_Report.txt');
        filename = ['ME_with_',num2str(params.N_me),'_time_steps.mat'];
        cd vary_time_steps
        save(filename,'params','me_results');
        cd ../
    end

end

if vary_rank_and_storage

    mkdir vary_rank_and_storage/

    params = Read_Inputs_from_Xml();
    me_results = Read_Results('Modified_Euler_Cost_Report.txt');
    filename = ['ME_with_rank_',num2str(params.rank),'_and_storage_',num2str(params.Maximum_Block_Update_Storage),'.mat'];
    cd vary_rank_and_storage/
    save(filename,'params','me_results');
    cd ../

end
