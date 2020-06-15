tic
model_tc_1_data = load([pwd '\DETECT_Versions\Efficient_DETECT\Efficient_Outputs\output_6hourly.mat']);
model_tc_2_data = load([pwd '\26_Jan_2020_17_38_cctype.mat']);


for time = 1:732,
    disp('processing Iteration: ')
    disp(time);
    for source = 1:4,
%     disp('processing Iteration')
%     disp(time)
    model_tc_1_data_iter = model_tc_1_data.out.CO2type(source,1:101,time);
    model_tc_2_data_iter = model_tc_2_data.x(source,1:101,time);
    scatter(model_tc_2_data_iter,model_tc_1_data_iter,2)
    hold on
    end;
end;
hold off
toc
