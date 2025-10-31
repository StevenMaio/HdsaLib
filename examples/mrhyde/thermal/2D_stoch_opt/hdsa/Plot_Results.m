classdef Plot_Results < handle

    properties
        m
        n
        ens_size
        coord_x
        coord_y
        time_steps
        X
        Y
    end

    methods
        function this = Plot_Results(ens_size)
            this.ens_size = ens_size;
            file_name = 'output_lofi_opt_sample_0.exo';
            this.coord_x = ncread(file_name, 'coordx');
            this.coord_y = ncread(file_name, 'coordy');
            this.time_steps = ncread(file_name, 'time_whole');
            this.m = length(this.time_steps);
            this.n = sqrt(size(this.coord_x,1));
            this.X = reshape(this.coord_x,this.n,this.n);
            this.Y = reshape(this.coord_y,this.n,this.n);
        end

        function [] = Plot_LoFi_Opt_State(this,ens_id)
            file_name = ['output_lofi_opt_sample_',num2str(ens_id),'.exo'];
            T = ncread(file_name,'vals_nod_var1');
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
    end
end