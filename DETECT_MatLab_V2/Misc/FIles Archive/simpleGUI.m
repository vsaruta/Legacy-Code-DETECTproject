function simpleGUI
    hFig = figure('Visible','off', 'Menu','none', 'Name','DETECT models', 'Resize','off', 'Position',[100 100 600 200]);
    movegui(hFig,'center')          %# Move the GUI to the center of the screen

    hBtnGrp = uibuttongroup('Position',[0 0 1 1], 'Units','Normalized');
    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[15 150 200 30], 'String','DETECT MatLab Original NSS', 'Tag','Model_1')
    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[15 120 200 30], 'String','DETECT MatLab Efficient NSS', 'Tag','Model_2')
    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[15  90 200 30], 'String','DETECT MatLab SS', 'Tag','Model_3')
    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[15  60 200 30], 'String','Divide', 'Tag','Model_4')

    uicontrol('Style','pushbutton', 'String','Run Model', 'Position',[15 30 60 25], 'Callback',{@button_callback})

%     hEdit1 = uicontrol('Style','edit', 'Position',[150 150 60 20], 'String','10');
%     hEdit2 = uicontrol('Style','edit', 'Position',[250 150 60 20], 'String','20');
    hEdit3 = uicontrol('Style','edit', 'Position',[100 30 60 25], 'String','');

    set(hFig, 'Visible','on')        %# Make the GUI visible

    %# callback function
    function button_callback(src,ev)
        switch get(get(hBtnGrp,'SelectedObject'),'Tag')
            case 'Model_1',
                res = 'Model_1';
%                 disp('works')
                Efficient_ScriptFile();
            case 'Model_2',  res = 'Model_2';
            case 'Model_3',  res = 'Model_3';
            case 'Model_4',  res = 'Model_4';
            otherwise, res = '';
        end
        set(hEdit3, 'String',res)
    end
end