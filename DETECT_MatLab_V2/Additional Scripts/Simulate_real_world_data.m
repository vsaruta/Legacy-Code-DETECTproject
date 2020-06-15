clear all
% simple_R_output = load('Outputs\Rout\cctype_imp_check5.mat');
simple_R_output = load([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly.mat'],'out')
time = 733;
source = 4;
depth = 101;
new_time = 1;
fake_data_time = NaN(source,depth,time);

% Here we set how long it take to get 1 messerment
% eg 24 hours - 1 meserment a day
% Min is 6 hours - 4 times a day
how_often = 24

if how_often < 6,
    how_often = 6;
end

n_messerments = how_often/6;
Total_N_messerments = 732/(how_often/6)


depth_i = 1:101;
for time_i = 1:733,
   for source_i = 1:4,
        while new_time ~= 733+n_messerments,
%             fake_data_time(source_i,depth_i,new_time) = simple_R_output.x(source_i,depth_i,new_time);
            fake_data_time(source_i,depth_i,new_time) = simple_R_output.out.CO2type(source_i,depth_i,new_time);
            new_time = new_time + n_messerments;
        end
    new_time = 1;
    end
end



fake_data_time_and_depth = NaN(source,depth,time);
depth_i = 1;
while depth_i ~= 111,
   fake_data_time_and_depth(:,depth_i,:) = fake_data_time(:,depth_i,:);
   depth_i = depth_i + 10;
end



% If needed we can change NaN cells to another type:
fake_data_time_and_depth_NaN = fake_data_time_and_depth;


% Change NaN to 0
fake_data_time_and_depth_0 = fake_data_time_and_depth_NaN;
fake_data_time_and_depth_0(isnan(fake_data_time_and_depth_0)) = 0;

% Change NaN to -9999
fake_data_time_and_depth_9999 = fake_data_time_and_depth_NaN;
fake_data_time_and_depth_9999(isnan(fake_data_time_and_depth_9999)) = -9999;

% Or fully remove NaN colums(we got 11 depth over 61 time steps for all 4 sources)
fake_data_time_and_depth_CLEAR = fake_data_time_and_depth;
fake_data_time_and_depth_CLEAR( :, all( isnan( fake_data_time_and_depth_CLEAR ), 1 ) ) = [];


% Or fully remove NaN colums(we got 11 depth over 61 time steps for all 4 sources)
% HERE NEED TO ADD ID BEFORE REMOVING
% fake_data_time_and_depth_CLEAR_ID = fake_data_time_and_depth;
% source
% time
% depth
% fake_data_time_and_depth_CLEAR_ID( :, all( isnan( fake_data_time_and_depth_CLEAR_ID ), 1 ) ) = [];

   
