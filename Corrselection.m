

function varargout = Corrselection(varargin)
% CORRSELECTION MATLAB code for Corrselection.fig
%      CORRSELECTION, by itself, creates a new CORRSELECTION or raises the existing
%      singleton*.
%
%      H = CORRSELECTION returns the handle to a new CORRSELECTION or the handle to
%      the existing singleton*.
%
%      CORRSELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRSELECTION.M with the given input arguments.
%
%      CORRSELECTION('Property','Value',...) creates a new CORRSELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Corrselection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Corrselection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Corrselection

% Last Modified by GUIDE v2.5 27-Jan-2016 09:57:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Corrselection_OpeningFcn, ...
                   'gui_OutputFcn',  @Corrselection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Corrselection is made visible.
function Corrselection_OpeningFcn(hObject, eventdata, handles, varargin)
global curveincl correlationcurves corfit meancorrcurve ItraceII correlationcurves2
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Corrselection (see VARARGIN)

% Choose default command line output for Corrselection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.slider3,'Max',size(correlationcurves,2)-1);
set(handles.slider3,'SliderStep',[1/(size(correlationcurves,2)-2) 1/(size(correlationcurves,2)-2)]);
correlationcurves2=correlationcurves;

for kl=2:size(correlationcurves,2)
if curveincl(kl-1)==0
    correlationcurves2(:,kl)=NaN;
end
end
meancorrcurve=nanmean(correlationcurves2(:,2:end),2);

set(handles.checkbox2,'Value',curveincl(get(handles.slider3,'Value'),1));
    plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
     hold (handles.axes2,'off')

     
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
hold (handles.axes1,'on')



semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
hold (handles.axes1,'off')

% UIWAIT makes Corrselection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Corrselection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
global curveincl correlationcurves corfit meancorrcurve ItraceII correlationcurves2
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


if get(hObject,'Value')<size(correlationcurves,2)
    
    correlationcurves2=correlationcurves;

for kl=2:size(correlationcurves,2)
if curveincl(kl-1)==0
    correlationcurves2(:,kl)=NaN;
end
end
meancorrcurve=nanmean(correlationcurves2(:,2:end),2);

set(handles.checkbox2,'Value',curveincl(get(hObject,'Value'),1));
get(hObject,'Min');
get(hObject,'Max');
     get(hObject,'Value')+1;
    plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(hObject,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(hObject,'Value')+1),1) mean(ItraceII(:,get(hObject,'Value')+1),1)],'r-')
     hold (handles.axes2,'off')

     
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(hObject,'Value')+1),'b+')
hold (handles.axes1,'on')



semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(hObject,'Value')+1),'k-','LineWidth',2)
axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(hObject,'Value')+1))) max(max(correlationcurves(1:10,get(hObject,'Value')+1)))]) 
hold (handles.axes1,'off')
end

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)

global curveincl correlationcurves correlationcurves2 meancorrcurve corfit ItraceII

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



curveincl(get(handles.slider3,'value'),1)=get(hObject,'Value');
correlationcurves2=correlationcurves;

for kl=2:size(correlationcurves,2)
if curveincl(kl-1)==0
    correlationcurves2(:,kl)=NaN;
end
end
meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
  
set(hObject,'Value',curveincl(get(handles.slider3,'Value'),1));
    plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
     hold (handles.axes2,'off')
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
hold (handles.axes1,'on')



semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
hold (handles.axes1,'off')


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global goon
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goon=1;
close Corrselection
