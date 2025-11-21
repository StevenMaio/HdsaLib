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

        function [] = Plot_Domain(this, save_filename)
            if nargin < 2
                save_filename = [];
            end

            % Create a new figure
            figure,

            % Set axis limits and aspect ratio
            axis([0 1 0 1]);
            axis equal;
            hold on;

            % Define the colors for the boundaries
            colors = lines(4);

            % Plot the boundaries with shading
            fill([0.0, 0.01, 0.01, 0.0], [.001, .001, .999, .999], colors(1,:), 'EdgeColor', 'none'); % Left boundary
            text(.03, 0.45, '$\partial \Omega_1$', 'Interpreter', 'latex', 'FontSize', 24);

            fill([.99, 1.0, 1.0, .99], [.001, .001, .999, .999], colors(2,:), 'EdgeColor', 'none'); % Right boundary
            text(0.85, 0.45, '$\partial \Omega_2$', 'Interpreter', 'latex', 'FontSize', 24);

            fill([.001, .99, .99, .001], [.99, .99, 1.0, 1.0], colors(3,:), 'EdgeColor', 'none'); % Top boundary
            text(0.45, 0.94, '$\partial \Omega_3$', 'Interpreter', 'latex', 'FontSize', 24);

            fill([.001, .99, .99, .001], [0.0, 0.0, 0.01, 0.01], colors(4,:), 'EdgeColor', 'none'); % Bottom boundary
            text(0.45, .05, '$\partial \Omega_4$', 'Interpreter', 'latex', 'FontSize', 24);

            % Create a shaded box in the region (0.4, 0.8)^2
            fill([0.4, 0.8, 0.8, 0.4], [0.4, 0.4, 0.8, 0.8], [0.5 0.5 0.5], 'EdgeColor', 'none'); % Gray box
            text(0.54, 0.6, '$f(\mathbf{\xi})$', 'Interpreter', 'latex', 'FontSize', 24);

            % Add an arrow pointing from left to right
            quiver(0.2, 0.2, 0.2, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 1, 'AutoScale', 'off');
            text(0.25, 0.2, 'v', 'FontSize', 24, 'VerticalAlignment', 'bottom');

            xticks([0,.2,.4,.6,.8,1])
            yticks([0,.2,.4,.6,.8,1])
            xlim([0,1])
            ylim([0,1])
            xlabel('x_1')
            ylabel('x_2')
            set(gca, 'FontSize', 24);
            if ~isempty(save_filename)
                print(save_filename, '-depsc');
            end
        end

        function [] = Plot_State_Comparison(this,num_ens, save_filename)
            if nargin < 3
                save_filename = [];
            end

            e = zeros(num_ens,11);
            for k = 1:num_ens
                d = this.Load_Discrepancy(k); 
                e(k,:) = vecnorm(d,2);
            end
            [~,i] = min(abs(e(:,end)-median(e(:,end))));

            [cmin_hifi,cmax_hifi] = this.Plot_HiFi_Opt_State(i,true);
            [cmin_lofi,cmax_lofi] = this.Plot_LoFi_Opt_State(i,true);

            cmin = min(cmin_lofi,cmin_hifi);
            cmax = max(cmax_lofi,cmax_hifi);

            clim([cmin,cmax])
            xlabel('x_1')
            ylabel('x_2')
            title('Low-fidelity State Solution $\tilde{S}(\tilde{z})$','Interpreter','latex')
            set(gca, 'FontSize', 24);
            pause(.1)
            if ~isempty(save_filename)
                print([save_filename,'_lofi'], '-depsc');
            end

            this.Plot_HiFi_Opt_State(i,true);
            clim([cmin,cmax])
            xlabel('x_1')
            ylabel('x_2')
            title('High-fidelity State Solution $S(\tilde{z})$','Interpreter','latex')
            set(gca, 'FontSize', 24);
            if ~isempty(save_filename)
                print([save_filename,'_hifi'], '-depsc');
            end

        end

        function [cmin,cmax] = Plot_State(this, T, final_time)
            if nargin < 3
                final_time = false;
            end

            if ~final_time

                cmin = min(T(:));
                cmax = max(T(:));

                figure,
                for k = 1:this.m
                    Z = reshape(T(:,k),this.n,this.n);
                    surf(this.X,this.Y,Z);
                    view(2)
                    title(['Time = ',num2str(this.time_steps(k))]);
                    clim([cmin,cmax]);
                    shading interp
                    colorbar();
                    %axis equal;
                    pause(.5);
                end

            else

                cmin = min(T(:,end));
                cmax = max(T(:,end));

                figure,
                Z = reshape(T(:,end),this.n,this.n);
                surf(this.X,this.Y,Z);
                view(2)
                title(['Time = ',num2str(this.time_steps(end))]);
                clim([cmin,cmax]);
                shading interp
                axis square
                xticks([0,.2,.4,.6,.8,1])
                yticks([0,.2,.4,.6,.8,1])
                xlim([0,1])
                ylim([0,1])
                colorbar();

            end
        end

        function [cmin,cmax] = Plot_HiFi_Opt_State(this,ens_id, final_time)
            if nargin < 3
                final_time = false;
            end
            T = this.Load_HiFi_State(ens_id);
            [cmin,cmax] = this.Plot_State(T, final_time);
        end

        function [cmin,cmax] = Plot_LoFi_Opt_State(this,ens_id, final_time)
            if nargin < 3
                final_time = false;
            end
            T = this.Load_LoFi_State(ens_id);
            [cmin,cmax] = this.Plot_State(T, final_time);
        end

        function [] = Plot_Discrepancy(this,ens_id)
            d = this.Load_Discrepancy(ens_id);
            this.Plot_State(d);
        end

        function [] = Plot_Optimal_z(this,num_samples, save_filename)
            if nargin < 3
                save_filename = [];
            end
            [z_mean, z_samples] = this.Load_Opt_z_Update(num_samples);
            z_lofi = load('optimization/final_params_.dat');
            z_hifi = load('hifi_optimization/final_params_hifi_.dat');

            s_mean = cell(4,1);
            s_lofi = cell(4,1);
            s_hifi = cell(4,1);
            for j = 1:4
                [x,s_mean{j}] = this.Map_Control_to_Mesh(z_mean,j);
                [~,s_lofi{j}] = this.Map_Control_to_Mesh(z_lofi,j);
                [~,s_hifi{j}] = this.Map_Control_to_Mesh(z_hifi,j);
            end

            if num_samples > 0

                s_samples = cell(4,1);
                for j = 1:4
                    sj = zeros(length(x),num_samples);
                    for k = 1:num_samples
                        [~,sj(:,k)] = this.Map_Control_to_Mesh(z_samples(:,k),j);
                    end
                    s_samples{j} = sj;
                end

                col = lines(3);

                figure,
                hold on
                plot(x,-s_lofi{1},'--','Color',col(1,:),'LineWidth',3)
                plot(x,-s_mean{1},'Color',col(2,:),'LineWidth',3)
                plot(x,-s_hifi{1},':','Color',col(3,:),'LineWidth',3)
                plot(x,-s_samples{1},'Color',[.9,.9,.9],'LineWidth',3)
                plot(x,-s_mean{1},'Color',col(2,:),'LineWidth',3)
                plot(x,-s_lofi{1},'--','Color',col(1,:),'LineWidth',3)
                plot(x,-s_hifi{1},':','Color',col(3,:),'LineWidth',3)
                legend({'LoFi','Mean','HiFi'},'Location','southwest')
                xlabel('y')
                ylabel('Controller')
                title('Left Boundary')
                set(gca, 'FontSize', 24);
                if ~isempty(save_filename)
                    print([save_filename,'_left'], '-depsc');
                end

                figure,
                hold on
                plot(x,-s_lofi{2},'--','Color',col(1,:),'LineWidth',3)
                plot(x,-s_mean{2},'Color',col(2,:),'LineWidth',3)
                plot(x,-s_hifi{2},':','Color',col(3,:),'LineWidth',3)
                plot(x,-s_samples{2},'Color',[.9,.9,.9],'LineWidth',3)
                plot(x,-s_mean{2},'Color',col(2,:),'LineWidth',3)
                plot(x,-s_lofi{2},'--','Color',col(1,:),'LineWidth',3)
                plot(x,-s_hifi{2},':','Color',col(3,:),'LineWidth',3)
                legend({'LoFi','Mean','HiFi'},'Location','southwest')
                xlabel('y')
                ylabel('Controller')
                title('Right Boundary')
                set(gca, 'FontSize', 24);
                if ~isempty(save_filename)
                    print([save_filename,'_right'], '-depsc');
                end


                figure,
                hold on
                plot(x,-s_lofi{3},'--','Color',col(1,:),'LineWidth',3)
                plot(x,-s_mean{3},'Color',col(2,:),'LineWidth',3)
                plot(x,-s_hifi{3},':','Color',col(3,:),'LineWidth',3)
                plot(x,-s_samples{3},'Color',[.9,.9,.9],'LineWidth',3)
                plot(x,-s_mean{3},'Color',col(2,:),'LineWidth',3)
                plot(x,-s_lofi{3},'--','Color',col(1,:),'LineWidth',3)
                plot(x,-s_hifi{3},':','Color',col(3,:),'LineWidth',3)
                legend({'LoFi','Mean','HiFi'},'Location','southwest')
                xlabel('x')
                ylabel('Controller')
                title('Top Boundary')
                set(gca, 'FontSize', 24);
                if ~isempty(save_filename)
                    print([save_filename,'_top'], '-depsc');
                end


                figure,
                hold on
                plot(x,-s_lofi{4},'--','Color',col(1,:),'LineWidth',3)
                plot(x,-s_mean{4},'Color',col(2,:),'LineWidth',3)
                plot(x,-s_hifi{4},':','Color',col(3,:),'LineWidth',3)
                plot(x,-s_samples{4},'Color',[.9,.9,.9],'LineWidth',3)
                plot(x,-s_mean{4},'Color',col(2,:),'LineWidth',3)
                plot(x,-s_lofi{4},'--','Color',col(1,:),'LineWidth',3)
                plot(x,-s_hifi{4},':','Color',col(3,:),'LineWidth',3)
                legend({'LoFi','Mean','HiFi'},'Location','southwest')
                xlabel('x')
                ylabel('Controller')
                title('Bottom Boundary')
                set(gca, 'FontSize', 24);
                if ~isempty(save_filename)
                    print([save_filename,'_bottom'], '-depsc');
                end


            else
                col = lines(3);

                figure,
                hold on
                plot(x,-s_lofi{1},'Color',col(1,:),'LineWidth',3)
                plot(x,-s_mean{1},'--','Color',col(2,:),'LineWidth',3)
                plot(x,-s_hifi{1},':','Color',col(3,:),'LineWidth',3)
                legend({'LoFi','Mean','HiFi'},'Location','southwest')
                xlabel('y')
                ylabel('Controller')
                title('Left Boundary')
                set(gca, 'FontSize', 24);
                if ~isempty(save_filename)
                    print([save_filename,'_left'], '-depsc');
                end

                figure,
                hold on
                plot(x,-s_lofi{2},'Color',col(1,:),'LineWidth',3)
                plot(x,-s_mean{2},'--','Color',col(2,:),'LineWidth',3)
                plot(x,-s_hifi{2},':','Color',col(3,:),'LineWidth',3)
                legend({'LoFi','Mean','HiFi'},'Location','southwest')
                xlabel('y')
                ylabel('Controller')
                title('Right Boundary')
                set(gca, 'FontSize', 24);
                if ~isempty(save_filename)
                    print([save_filename,'_right'], '-depsc');
                end

                figure,
                hold on
                plot(x,-s_lofi{3},'Color',col(1,:),'LineWidth',3)
                plot(x,-s_mean{3},'--','Color',col(2,:),'LineWidth',3)
                plot(x,-s_hifi{3},':','Color',col(3,:),'LineWidth',3)
                legend({'LoFi','Mean','HiFi'},'Location','southwest')
                xlabel('x')
                ylabel('Controller')
                title('Top Boundary')
                set(gca, 'FontSize', 24);
                if ~isempty(save_filename)
                    print([save_filename,'_top'], '-depsc');
                end

                 figure,
                hold on
                plot(x,-s_lofi{3},'Color',col(1,:),'LineWidth',3)
                plot(x,-s_mean{3},'--','Color',col(2,:),'LineWidth',3)
                plot(x,-s_hifi{3},':','Color',col(3,:),'LineWidth',3)
                legend({'LoFi','Mean','HiFi'},'Location','southwest')
                xlabel('x')
                ylabel('Controller')
                title('Bottom Boundary')
                set(gca, 'FontSize', 24);
                if ~isempty(save_filename)
                    print([save_filename,'_bottom'], '-depsc');
                end
            end
        end

        function [x,s] = Map_Control_to_Mesh(this,z, component)
            
            n_grid = 100;
            x = linspace(.2,.8,n_grid)';
            offset = 5*(component-1);
            
            s = zeros(n_grid,1);
            for k = 1:5
                c = .3 + (k-1)*.1;
                s = s + z(offset+k) * exp(-400*(x-c).^2);
            end
        end

        function [] = Plot_Ensemble_Set(this, save_filename)
            if nargin < 2
                save_filename = [];
            end

            [sample_set, sample_weights] = this.Load_Sample_Data();
            ens_size = size(sample_set,1);
            
            figure,
            hold on
            I = find(sample_weights<.6);
            for k = 1:length(I)
                plot(sample_set(I(k),1), sample_set(I(k),2),'x','MarkerSize',15,'Color','black')
            end
            I = find(sample_weights>.6);
            for k = 1:length(I)
                plot(sample_set(I(k),1), sample_set(I(k),2),'d','MarkerSize',15,'Color','red')
            end
            xlim([.3,.9])
            ylim([.3,.9])
            xticks([0,.2,.4,.6,.8,1])
            yticks([0,.2,.4,.6,.8,1])
            xlabel('x_1')
            ylabel('x_2')
            title('Low-fidelity')
            axis square
            set(gca, 'FontSize', 24);
            if ~isempty(save_filename)
                print(save_filename, '-depsc');
            end

            [sample_set, sample_weights] = this.Load_Sample_Data_HiFi();
            figure,
            hold on
            I = find(sample_weights<.6);
            for k = 1:length(I)
                plot(sample_set(I(k),1), sample_set(I(k),2),'x','MarkerSize',15,'Color','black')
            end
            I = find(sample_weights>.6);
            for k = 1:length(I)
                plot(sample_set(I(k),1), sample_set(I(k),2),'d','MarkerSize',15,'Color','red')
            end
            xlim([.3,.9])
            ylim([.3,.9])
            xticks([0,.2,.4,.6,.8,1])
            yticks([0,.2,.4,.6,.8,1])
            xlabel('x_1')
            ylabel('x_2')
            title('High-fidelity')
            axis square
            set(gca, 'FontSize', 24);
            if ~isempty(save_filename)
                print([save_filename,'_hifi'], '-depsc');
            end

        end

        function [] = Plot_Prior_Delta_z_opt(this, ens_id, sample_id, final_time)
            if nargin < 4
                final_time = false;
            end
            prior_delta_z_opt = this.Load_Prior_Delta_z_opt(ens_id, sample_id);
            this.Plot_State(prior_delta_z_opt, final_time);
        end

        function [] = Plot_Time_Evolution(this)
            prior_delta_z_opt_time_evol_1 = load('hdsa/hdsa_output/prior/summary_statistics/prior_delta_z_opt_time_evol_1.txt');
            prior_discrep_data_time_evol_1 = load('hdsa/hdsa_output/prior/summary_statistics/prior_discrep_data_time_evol_1.txt');
            figure,
            hold on
            plot(this.time_steps,prior_delta_z_opt_time_evol_1,'Color',.9*ones(3,1))
            plot(this.time_steps,prior_discrep_data_time_evol_1,'Color','red')
        end

        function [i_min, j_min, i_max, j_max, D, sample_set] = Plot_Sample_Distances(this)
            sample_set = this.Load_Sample_Data();
            ens_size = size(sample_set,1);
            D = zeros(ens_size,ens_size);
            for i = 1:ens_size
                for j = 1:ens_size
                    D(i,j) = norm(sample_set(i,:)-sample_set(j,:));
                end
            end

            A = D + eye(ens_size);
            [~, linearIndex] = min(A(:));
            [i_min, j_min] = ind2sub(size(A), linearIndex);

            [~, linearIndex] = max(D(:));
            [i_max, j_max] = ind2sub(size(A), linearIndex);
        end

        function [sample_set, sample_weights] = Load_Sample_Data(this)
            sample_set = load('optimization/sample_set.dat');
            sample_weights = load('optimization/sample_weights.dat');
        end

        function [sample_set, sample_weights] = Load_Sample_Data_HiFi(this)
            sample_set = load('hifi_optimization/sample_set.dat');
            sample_weights = load('hifi_optimization/sample_weights.dat');
        end

        function [prior_delta_z_opt] = Load_Prior_Delta_z_opt(this, ens_id, sample_id)
            file_name = ['hdsa/hdsa_output/prior/prior_delta_z_opt/Vector_',num2str(sample_id),'_ens_',num2str(ens_id),'.exo'];
            prior_delta_z_opt = ncread(file_name,'vals_nod_var1');
        end

        function [prior_delta_z_pert] = Load_Prior_Delta_z_pert(this, ens_id, sample_id)
            file_name = ['hdsa/hdsa_output/prior/prior_delta_z_pert_1/Vector_',num2str(sample_id),'_ens_',num2str(ens_id),'.exo'];
            prior_delta_z_pert = ncread(file_name,'vals_nod_var1');
        end

        function [d] = Load_Discrepancy(this,ens_id)
            d = this.Load_HiFi_State(ens_id) - this.Load_LoFi_State(ens_id);
        end

        function [T] = Load_HiFi_State(this,ens_id)
            file_name = ['fwd_hifi_data_generation/output_hifi_sample_',num2str(ens_id-1),'.exo'];
            T = ncread(file_name,'vals_nod_var1');
        end

        function [T] = Load_LoFi_State(this,ens_id)
            file_name = ['optimization/output_lofi_opt_sample_',num2str(ens_id-1),'.exo'];
            T = ncread(file_name,'vals_nod_var1');
        end

        function [z_mean, z_samples] = Load_Opt_z_Update(this,num_samples)
            z_samples = zeros(20,num_samples);
            for k = 1:num_samples
                file_name = ['hdsa/hdsa_output/posterior/z_update/posterior_samples/Vector_',num2str(k),'.txt'];
                z_samples(:,k) = load(file_name);
            end
            z_mean = load('hdsa/hdsa_output/posterior/z_update/mean.txt');
        end

    end
end