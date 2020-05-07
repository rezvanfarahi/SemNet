function createfigure(cdata1)
%CREATEFIGURE(CDATA1)
%  CDATA1:  image cdata

%  Auto-generated by MATLAB on 28-Nov-2017 14:26:58

% Create figure
figure1 = figure;
colormap('hot');

% Create axes
axes1 = axes('Parent',figure1,'Layer','top','YDir','reverse');
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.5 74.5]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.5 74.5]);
box(axes1,'on');
hold(axes1,'on');

% Create image
image(cdata1,'Parent',axes1,'CDataMapping','scaled');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.265134275618375 0.694915254237288 0.0846890459363958 0.0866290018832383],...
    'Color',[1 1 1],...
    'LineWidth',2,...
    'LineStyle',':');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.349056537102473 0.604519774011299 0.0855724381625442 0.0903954802259874],...
    'Color',[1 1 1],...
    'LineWidth',2,...
    'LineStyle',':');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.433862190812721 0.536723163841808 0.0626042402826856 0.0677966101694908],...
    'Color',[1 1 1],...
    'LineWidth',2,...
    'LineStyle',':');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.49569964664311 0.406779661016949 0.126208480565371 0.131826741996233],...
    'Color',[1 1 1],...
    'LineWidth',2,...
    'LineStyle',':');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.621141342756184 0.382297551789077 0.0219681978798587 0.0263653483992467],...
    'Color',[1 1 1],...
    'LineWidth',2,...
    'LineStyle',':');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.643226148409894 0.210922787193974 0.167727915194346 0.175141242937853],...
    'Color',[1 1 1],...
    'LineWidth',2,...
    'LineStyle',':');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.810187279151943 0.109227871939736 0.0944063604240283 0.0998116760828624],...
    'Color',[1 1 1],...
    'LineWidth',2,...
    'LineStyle',':');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.129975265017668 0.779661016949152 0.135925795053004 0.145009416195857],...
    'Color',[1 1 1],...
    'LineWidth',2,...
    'LineStyle',':');

