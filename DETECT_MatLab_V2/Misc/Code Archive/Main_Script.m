finished = false;
while finished ~= true,
    fprintf('\n \n \nWelcome to DETECT models initialization script. \n')
    modeling_or_ploting = false;
    while modeling_or_ploting ~= true,
       what_to_do = input('You can choose to: \n 1 - Run the models \n 2 - Produce Plots \n 3 - Quit \n Whad do you choose to do: \n');
       if what_to_do == 1,
           fprintf('\n\nCurrently such model(s) output(s) already exist(s) in the system: \n')
           count = 0;
           if isfile([strrep(pwd,'\','/') '/DETECT_Versions/Original_DETECT/Outputs/output_6hourly'  '.mat']),
               count = count + 1;
               fprintf('Original_DETECT\n')
           end
           
           if isfile([strrep(pwd,'\','/') '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly'  '.mat']),
               count = count + 1;
               fprintf('Efficient_DETECT \n')
           end
           if count == 0,
               fprintf('None \n')
           end
           
           selected = false;
           while selected ~= true,
               model_to_Run = input('You can choose to run: \n 1 - Original DETECT \n 2 - Efficient DETECT \n 3 - Cancel \n Which model to run: \n');
               if model_to_Run == 1,
                   ScriptFile()
                   selected = true;
               elseif model_to_Run == 2,
                   Efficient_ScriptFile()
                   selected = true;
               elseif model_to_Run == 3,
                   selected = true;
               else
                   fprintf('ERROR: No such model. Please, choose correct number! \n');
               end
           end
           modeling_or_ploting = true;
           
       elseif what_to_do == 2,
            correct_select = false;
            while correct_select ~= true,
                to_plot_or_not = input('Would you like to plot model output or compare models? \n 1 - Plot Output \n 2 - Compare Models \n Choose 1 or 2: \n');
                if to_plot_or_not == 1,
                    fprintf('\n\nCurrently such model(s) output(s) already exist(s) in the system: \n')
                    count = 0;
                    if isfile([strrep(pwd,'\','/') '/DETECT_Versions/Original_DETECT/Outputs/output_6hourly'  '.mat']),
                        count = count + 1;
                        fprintf('Original_DETECT\n')
                    end

                    if isfile([strrep(pwd,'\','/') '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly'  '.mat']),
                        count = count + 1;
                        fprintf('Efficient_DETECT \n')
                    end
                    if count == 0,
                        fprintf('None \n')
                    end
%                     HERE NEED TO FINISH DEVELOPMENT
                    result_of_model = input('You can choose to plot results of: \n 1 - Original DETECT \n 2 - Efficient DETECT \n Which model output to use: \n');
                    if result_of_model == 1,
                        current_plot_data = load([pwd '/DETECT_Versions/Original_DETECT/Outputs/output_6hourly'  '.mat'],'out')
                        plot(current_plot_data.out.CO2)     
                    elseif result_of_model == 2,
%                         current_plot_data = load([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly'  '.mat'],'out')
%                         plot(current_plot_data)
                        Efficient_Rsoil_timeseries()
                    else
                       fprintf('ERROR: Sorry! There is no DATA to plot. Plotting Original DETECT \n')
                       current_plot_data = load([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly'  '.mat'],'out')
                    end
%                     plot(current_plot_data.out.CO2)
                    correct_select = true;
                elseif to_plot_or_not == 2,
                    model_array_1 = load([pwd '/DETECT_Versions/Original_DETECT/Outputs/output_6hourly'  '.mat'],'out')
                    model_array_2 = load([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly'  '.mat'],'out')
                    
                    to_compare_1 = input('You can choose to run: \n 1 - Original DETECT \n 2 - Efficient DETECT \n 3 - Cancel \n Which model to run: \n');
                    to_compare_2 = input('You can choose to run: \n 1 - Original DETECT \n 2 - Efficient DETECT \n 3 - Cancel \n Which model to run: \n');
                    if to_compare_1 == 1,
                        first_model = model_array_1.out.CO2;
                    elseif to_compare_1 ==2,
                        first_model = model_array_2.out.CO2;
                    end
                    if to_compare_1 == 1,
                        second_model = model_array_1.out.CO2;
                    elseif to_compare_1 ==2,
                        second_model = model_array_2.out.CO2;
                    end
%                     scatterplot(float(first_model(1:100,1)),float(second_model(1:100,1)))
                    fprintf('Sorry! Currently unsuported.')
                    correct_select = true;
                else
                    fprintf('ERROR: No such option. Please, choose correct number "1" or "2"! \n')
                end

                % DETECT_Efficient_output = load([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly'  '.mat'],'out')
                % plot(DETECT_Efficient_output)
            end
           modeling_or_ploting = true;
           
       elseif what_to_do == 3,
            finished = true;
            modeling_or_ploting = true;
       else
           fprintf('ERROR: No such option. Please, choose correct answear "1","2", or "3"! \n')
       end 
    end
    

    correct_finished = false;
    while correct_finished ~= true,
        continue_or_finish = input('\n\nDo you want to do something else or finish? \n 1 - Continue \n 2 - Finish \n Choose 1 or 2: \n');
        if continue_or_finish == 1,
            correct_finished = true;
            finished = false;
        elseif continue_or_finish == 2,
            correct_finished = true;
            finished = true;
             fprintf('\n\n\nThank you for using DETECT initialization script. \nHave a nice day! \n')
        else
             fprintf('ERROR: No such option. Please, choose correct number "1" or "2"! \n')
        end
    end
end