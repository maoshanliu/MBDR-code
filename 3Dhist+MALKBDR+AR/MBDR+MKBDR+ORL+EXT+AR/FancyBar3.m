function [ax] = FancyBar3(varargin)
% 3D histogram
% author : xiaoyeyimier
% date   : 2022-3-30
% version : 1.01
% contact : xiaoyeyimier@163.com
%=========================================================================%
%% Input
% X          : The abscissa of each set of data.
%              Each row of data must be uniform and increasing
% Y          : The position of each histogram on the plane
% Z          : yData of the data of each group
% B3settings : Drawing parameters,
% |            including the length, width and height of each cube
% |___dWidth         : The interval between adjacent pieces of data in the 
% |                    same set;
% |___LineWidth      : The border thickness of the histogram;
% |___TempYWidth     : The interval between different datasets;
% |___TempFaceAlpha  : The opacity of the histogram;
% |___ViewAngle      : angle of view ;
% |___colorarr       : Colors used for different groups of data;
% ax          : You can specify a coordinate handle that is currently in use
%=========================================================================%
%% Output
% ax         : Handle to the coordinate axes used for drawing
%=========================================================================%
% read & check parameters
[X,Y,Z,B3,ax] = inputParamCheck(varargin,nargin);
B3 = StructCompletion(B3,X);
TempFaceAlpha  = B3.TempFaceAlpha;
colorarr = B3.colorarr;
LineWidth = B3.LineWidth;
% draw plane by plane
for looptimeI = 1:length(Y)
    if looptimeI==length(Y)
        TempYWidth = B3.TempYWidth*(Y(looptimeI)-Y(looptimeI-1));
    else
        TempYWidth = B3.TempYWidth*(Y(looptimeI+1)-Y(looptimeI));
    end
    Xdata = X(looptimeI,:);
    Zdata = Z(looptimeI,:);
    Ydata = Y(looptimeI);
    TempWidth = (X(looptimeI,2)-X(looptimeI,1))*B3.dWidth;
    for looptimeJ = 1:length(Xdata)
        TempXdata = Xdata(looptimeJ);
        TempZdata = Zdata(looptimeJ);
        fillres = fill3([TempXdata-TempWidth/2,TempXdata-TempWidth/2,...
            TempXdata+TempWidth/2,TempXdata+TempWidth/2],...
            [Ydata-TempYWidth/2,Ydata-TempYWidth/2,Ydata-TempYWidth/2,Ydata-TempYWidth/2],...
            [0,TempZdata,TempZdata,0],colorarr(looptimeI,:));
        fillres.FaceAlpha = TempFaceAlpha;
        fillres.LineWidth = LineWidth;
        fillres = fill3([TempXdata-TempWidth/2,TempXdata-TempWidth/2,...
            TempXdata+TempWidth/2,TempXdata+TempWidth/2],...
            [Ydata+TempYWidth/2,Ydata+TempYWidth/2,Ydata+TempYWidth/2,Ydata+TempYWidth/2],...
            [0,TempZdata,TempZdata,0],colorarr(looptimeI,:));
        fillres.FaceAlpha = TempFaceAlpha;
        fillres.LineWidth = B3.LineWidth;
        fillres = fill3([TempXdata-TempWidth/2,TempXdata-TempWidth/2,...
            TempXdata-TempWidth/2,TempXdata-TempWidth/2],...
            [Ydata+TempYWidth/2,Ydata+TempYWidth/2,Ydata-TempYWidth/2,Ydata-TempYWidth/2],...
            [0,TempZdata,TempZdata,0],colorarr(looptimeI,:));
        fillres.FaceAlpha = TempFaceAlpha;
        fillres.LineWidth = B3.LineWidth;
        fillres = fill3([TempXdata+TempWidth/2,TempXdata+TempWidth/2,...
            TempXdata+TempWidth/2,TempXdata+TempWidth/2],...
            [Ydata+TempYWidth/2,Ydata+TempYWidth/2,Ydata-TempYWidth/2,Ydata-TempYWidth/2],...
            [0,TempZdata,TempZdata,0],colorarr(looptimeI,:));
        fillres.FaceAlpha = TempFaceAlpha;
        fillres.LineWidth = B3.LineWidth;
        fillres = fill3([TempXdata-TempWidth/2,TempXdata+TempWidth/2,...
            TempXdata+TempWidth/2,TempXdata-TempWidth/2],...
            [Ydata-TempYWidth/2,Ydata-TempYWidth/2,Ydata+TempYWidth/2,Ydata+TempYWidth/2],...
            [TempZdata,TempZdata,TempZdata,TempZdata],colorarr(looptimeI,:));
        fillres.FaceAlpha = TempFaceAlpha;
        fillres.LineWidth = B3.LineWidth;
        hold on;
    end
end
view(B3.ViewAngle);
xlim(ax,[min(min(X)),max(max(X))]);
grid on;
deltaX = (max(max(X))-min(min(X)))*0.05;
xlim([min(min(X))-TempWidth/2-deltaX,max(max(X))+TempWidth/2+deltaX]);
deltaY = (max(max(Y))-min(min(Y)))*0.05;
ylim([min(min(Y))-TempYWidth/2-deltaY,max(max(Y))+TempYWidth/2+deltaY]);
end

function [X,Y,Z,B3settings,ax] = inputParamCheck(V,N)
if N==5
    X = V{1};
    Y = V{2};
    Z = V{3};
    ax = V{4};
    B3settings = V{5};
elseif N==4
    X = V{1};
    Y = V{2};
    Z = V{3};
    ax = V{4};
    B3settings = [];
elseif N==3
    X = V{1};
    Y = V{2};
    Z = V{3};
    ax = gca;
    B3settings = [];
else
    error('Check the number of input parameters')
end
end

function B3C = StructCompletion(B3,X)
if ~isstruct(B3)
    B3C.dWidth         = 0.4;
    B3C.LineWidth      = 0.7;
    B3C.TempYWidth     = 0.3;
    B3C.TempFaceAlpha  = 0.7;
    B3C.ViewAngle      = [45,45];
    B3C.colorarr       = winter(size(X,1));
else
    B3C = B3;
    if ~isfield(B3C,'dWidth')
        B3C.dWidth = 0.4;
    end
    if ~isfield(B3C,'LineWidth')
        B3C.LineWidth = 0.7;
    end
    if ~isfield(B3C,'TempYWidth')
        B3C.TempYWidth = 0.3;
    end
    if ~isfield(B3C,'TempFaceAlpha')
        B3C.TempFaceAlpha = 0.7;
    end
    if ~isfield(B3C,'colorarr')
        B3C.colorarr = winter(size(X,1));
    end
    if ~isfield(B3C,'ViewAngle')
        B3C.ViewAngle = [45,45];
    end
end
end