% fig = uifigure('Position',[680 678 398 271]);
% bg = uibuttongroup(fig,'Position',[137 113 123 85]); 
% btn = uibutton(fig);
% rb1 = uiradiobutton(bg,'Position',[10 60 91 15]);
% rb2 = uiradiobutton(bg,'Position',[10 38 91 15]);
% rb3 = uiradiobutton(bg,'Position',[10 16 91 15]);
% 
% 
% rb1.Text = 'Model 1';
% rb2.Text = 'Model 2';
% rb3.Text = 'Model 3';
% 
% 
% % rb3.Value = true;
% % fig = uifigure;
% 
% 
% 
% if bg.rb1.Value == 1
%   % code
%   dispt('N1 selected')
% elseif bg.rb2.Value == 1
%   % other code
%   disp('N2 selected')
% else
%     disp('oops')
% end


function buttonPlot
% Create a figure window
fig = uifigure;

% Create a UI axes
ax = uiaxes('Parent',fig,...
            'Units','pixels',...
            'Position', [104, 123, 300, 201]);   

% Create a push button
btn = uibutton(fig,'push',...
               'Position',[420, 218, 100, 22],...
               'ButtonPushedFcn', @(btn,event) plotButtonPushed(btn,ax));
end

% Create the function for the ButtonPushedFcn callback
function plotButtonPushed(btn,ax)
        x = linspace(0,2*pi,100);
        y = sin(x);
        plot(ax,x,y)
end