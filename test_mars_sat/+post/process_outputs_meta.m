function sim_process = process_outputs_meta(sim)

sim_process = sim;

fail_percent = 0;
id_time = 0;
id_time_count = 0;
fail_list = [];
x1_all = [];
x2_all = [];
unknown_list = [];
for i = 1:sim.config.N_trials
    data = sim.outputs{i};
    fail_percent = fail_percent + data.fail;
    if(data.fail == 1)
        fail_list(end+1) = i;
    end
    if(data.true_mode == -1)
        unknown_list(end+1) = i;
    end
    if(data.true_mode >= 0 && data.true_mode == data.id_mode)
        id_time = id_time + data.k_end;
        id_time_count = id_time_count + 1;
    end
    x1_all = [x1_all, data.x_hist(1,:)];
    x2_all = [x2_all, data.x_hist(2,:)];
end
fail_percent = 100*(fail_percent/sim.config.N_trials);

sim_process.outputs_meta.fail_percent = fail_percent;
sim_process.outputs_meta.success_id_time_for_known_fault = id_time/id_time_count; 
sim_process.outputs_meta.fail_list = fail_list;
sim_process.outputs_meta.unknown_list = unknown_list;
sim_process.outputs_meta.x1_all = x1_all;
sim_process.outputs_meta.x2_all = x2_all;
% sim_process.outputs_meta.bool_violate_constraints = check_constraints(fid, x1_all, x2_all);
% sim_process.outputs_meta.violate_constraints_percent = 100.0*mean(sim_process.outputs_meta.bool_violate_constraints, 2);
% sim_process.outputs_meta.violate_percent_total = sum(sim_process.outputs_meta.violate_constraints_percent);

end

% function bool_violate_constraints = check_constraints(fid, x1_all, x2_all)
% bool_violate_constraints = zeros(3, length(x1_all));
% for i = 1:length(x1_all)
%     if(x1_all(i) > fid.config.x1_max)
%         bool_violate_constraints(1,i) = 1;
%         break;
%     end
%     if(x1_all(i) < x2_all(i))
%         bool_violate_constraints(2,i) = 1;
%         break;
%     end
%     if(x2_all(i) < fid.config.x2_min)
%         bool_violate_constraints(3,i) = 1;
%         break;
%     end
% end
% end