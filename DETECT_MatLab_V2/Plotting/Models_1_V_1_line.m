tic
model_tc_1_data = load([pwd '\DETECT_Versions\Original_DETECT\Outputs\output_6hourly.mat']);
model_tc_2_data = load([pwd '\DETECT_Versions\Efficient_DETECT\Efficient_Outputs\output_6hourly.mat']);


for time = 1:732,
    for source = 1:4,
    model_tc_1_data_iter = model_tc_1_data.out.CO2type(source,1:101,time);
    model_tc_2_data_iter = model_tc_2_data.out.CO2type(source,1:101,time);
    scatter(model_tc_2_data_iter,model_tc_1_data_iter,2)
    hold on
    end;
end;
hold off
toc
