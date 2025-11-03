classdef Plot_Results < handle

    properties
        m
        n 
        coord_x
        coord_y
        time_steps
        X
        Y
    end

    methods
        function this = Plot_Results()
            file_name = 'optimization/output_lofi_opt_sample_0.exo';
            this.coord_x = ncread(file_name, 'coordx');
            this.coord_y = ncread(file_name, 'coordy');
            this.time_steps = ncread(file_name, 'time_whole');
            this.m = length(this.time_steps);
            this.n = sqrt(size(this.coord_x,1));
            this.X = reshape(this.coord_x,this.n,this.n);
            this.Y = reshape(this.coord_y,this.n,this.n);
        end

        function [] = Plot_State(this,T)
            cmin = min(T(:));
            cmax = max(T(:));
            close all
            figure(1)
            for k = 1:this.m
                Z = reshape(T(:,k),this.n,this.n);
                figure(1)
                surf(this.X,this.Y,Z);
                view(2)
                title(['Time = ',num2str(this.time_steps(k))]);
                clim([cmin,cmax]);
                colorbar();
                pause(.5);
            end
        end

        function [] = Plot_LoFi_Opt_State(this,ens_id)
            T = this.Load_LoFi_State(ens_id);
            this.Plot_State(T);
        end

        function [] = Plot_Discrepancy(this,ens_id)
            d = this.Load_Discrepancy(ens_id);
            this.Plot_State(d);
        end

        function [] = Plot_Optimal_z(this,num_samples)
            [z_mean, z_samples] = this.Load_Opt_z_Update(num_samples);
            z_lofi = load('optimization/final_params_.dat');
            z_hifi = load('hifi_optimization/final_params_hifi_.dat');

            figure,
            hold on
            for k = 1:num_samples
                plot(1:4,-z_samples(:,k),'o','MarkerSize',10,'Color','black')
            end
            plot(1:4,-z_mean,'x','MarkerSize',10,'Color','black')
            plot(1:4,-z_lofi,'s','MarkerSize',10,'Color','red')
            plot(1:4,-z_hifi,'d','MarkerSize',10,'Color','magenta')
            xlim([0,5])
        end

        function [] = Plot_Ensemble_Set(this)
            [sample_set, sample_weights] = this.Load_Sample_Data();
            ens_size = size(sample_set,1);
            figure,
            hold on
            for k = 1:ens_size
                plot(sample_set(k,1), sample_set(k,2),'x','MarkerSize',10,'Color','black')
            end
            I = find(sample_weights>.6);
            for k = 1:length(I)
                plot(sample_set(I(k),1), sample_set(I(k),2),'x','MarkerSize',10,'Color','red')
            end
            xlim([.3,.9])
            ylim([.3,.9])
        end

        function [sample_set, sample_weights] = Load_Sample_Data(this)
            sample_set = load('optimization/sample_set.dat');
            sample_weights = load('optimization/sample_weights.dat');
        end

        function [d] = Load_Discrepancy(this,ens_id)
            d = this.Load_HiFi_State(ens_id) - this.Load_LoFi_State(ens_id);
        end

        function [T] = Load_HiFi_State(this,ens_id)
            file_name = ['fwd_hifi_data_generation/output_hifi_sample_',num2str(ens_id),'.exo'];
            T = ncread(file_name,'vals_nod_var1');
        end

        function [T] = Load_LoFi_State(this,ens_id)
            file_name = ['optimization/output_lofi_opt_sample_',num2str(ens_id),'.exo'];
            T = ncread(file_name,'vals_nod_var1');
        end

        function [z_mean, z_samples] = Load_Opt_z_Update(this,num_samples)
            z_samples = zeros(4,num_samples);
            for k = 1:num_samples
                file_name = ['hdsa/hdsa_output/posterior/z_update/posterior_samples/Vector_',num2str(k),'.txt'];
                z_samples(:,k) = load(file_name);
            end
            z_mean = load('hdsa/hdsa_output/posterior/z_update/mean.txt');
        end

    end
end