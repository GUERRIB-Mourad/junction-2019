function varargout = testbigtp(varargin)
% TESTBIGTP M-file for testbigtp.fig
%      TESTBIGTP, by itself, creates a new TESTBIGTP or raises the existing
%      singleton*.
%
%      H = TESTBIGTP returns the handle to a new TESTBIGTP or the handle to
%      the existing singleton*.
%
%      TESTBIGTP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTBIGTP.M with the given input arguments.
%
%      TESTBIGTP('Property','Value',...) creates a new TESTBIGTP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testbigtp_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testbigtp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testbigtp

% Last Modified by GUIDE v2.5 16-Nov-2019 22:19:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testbigtp_OpeningFcn, ...
                   'gui_OutputFcn',  @testbigtp_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before testbigtp is made visible.
function testbigtp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testbigtp (see VARARGIN)

% Choose default command line output for testbigtp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using testbigtp.
if strcmp(get(hObject,'Visible'),'off')
    plot(cos(0));
end

% UIWAIT makes testbigtp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = testbigtp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in tracer.
function tracer_Callback(hObject, eventdata, handles)
% hObject    handle to tracer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1 
        Ndebut=10;    %enlever les 10 premier echantillons
[fname,fpath]=uigetfile('*.wav');    %selection du fichier audio a ecouter
[x, fe, nbits]=wavread([fpath,fname]);   %chargement du fichier echantillonées en y avec la frequence et le nombre de bits
x=x/max(abs(x));     %...normalise
sound(x,fe);    %ecouter les echantillons avec la frequence donnée
s=length(x);    %longueur du son...
X = dct(x);
[XX,ind] = sort(abs(X),'descend');
need = 1;
s= eval(get(handles.edit1,'string'));
while norm(X(ind(1:need)))/norm(X)<s
   need = need+1;
end

xpc = need/length(X)*100;
X(ind(need+1:end)) = 0;
xx = idct(X);
%R=[x.-xx]
%plot([x;xx;x-xx]')
plot(x,'b');
hold on
plot(xx,'r')
hold on
plot(x-xx,'y')
legend('Original',[int2str(xpc) '% de coeffs.'],'Différence','Location','best') %....



    case 2
      % Load two speech waveforms of the same utterance (from TIMIT)
[fname,fpath]=uigetfile('*.wav');    %selection du fichier audio a ecouter
[d1,sr]=wavread([fpath,fname]);  
[fname,fpath]=uigetfile('*.wav');    %selection du fichier audio a ecouter
[d2,sr]=wavread([fpath,fname]);  
 
 % Listen to them together:
 ml = min(length(d1),length(d2));
 soundsc(d1(1:ml)+d2(1:ml),sr)
 % or, in stereo
 soundsc([d1(1:ml),d2(1:ml)],sr)

 % Calculate STFT features for both sounds (25% window overlap)
 D1 = specgram(d1,512,sr,512,384);
 D2 = specgram(d2,512,sr,512,384);

 % Construct the 'local match' scores matrix as the cosine distance 
 % between the STFT magnitudes
 ED1 = sqrt(sum((abs(D1)).^2));
 ED2 = sqrt(sum((abs(D2)).^2));
 SM = ((abs(D1))'*(abs(D2)))./(ED1'*ED2);
%%%%%%%%%%  SM = simmx(abs(D1),abs(D2));
 % Look at it:
%  subplot(121)
%  imagesc(SM)
  colormap(1-gray)
 % You can see a dark stripe (high similarity values) approximately
 % down the leading diagonal.

 % Use dynamic programming to find the lowest-cost path between the 
 % opposite corners of the cost matrix
 % Note that we use 1-SM because dp will find the *lowest* total cost
 M=1-SM;
 [r,D2x] = size(M);

% costs
C = zeros(r+1, D2x+1);
C(1,:) = NaN;
C(:,1) = NaN;
C(1,1) = 0;
C(2:(r+1), 2:(D2x+1)) = M;

% traceback
phi = zeros(r,D2x);

for i = 1:r; 
  for j = 1:D2x;
    [dmax, tb] = min([C(i, j), C(i, j+1), C(i+1, j)]);
    C(i+1,j+1) = C(i+1,j+1)+dmax;
    phi(i,j) = tb;
  end
end

% Traceback from top left
i = r; 
j = D2x;
p = i;
q = j;
while i > 1 & j > 1
  tb = phi(i,j);
  if (tb == 1)
    i = i-1;
    j = j-1;
  elseif (tb == 2)
    i = i-1;
  elseif (tb == 3)
    j = j-1;
  else    
    error;
  end
  p = [i,p];
  q = [j,q];
end

% Strip off the edges of the D matrix before returning
C = C(2:(r+1),2:(D2x+1));
 %%%%%%%%%%[p,q,C] = dp(1-SM);
 % Overlay the path on the local similarity matrix
 hold on; plot(q,p,'r'); hold off
 % Path visibly follows the dark stripe
 
 % Plot the minimum-cost-to-this point matrix too
 %subplot(122)
 imagesc(C)
 hold on; plot(q,p,'r'); hold off
%[image of DTW path]
 
 % Bottom right corner of C gives cost of minimum-cost alignment of the two
 M=C(size(C,1),size(C,2))
 
 if M<2
     axes(handles.axes2);
     cla;
     plot(D1);
     axes(handles.axes3);
     cla;
     plot(M);
 else
     axes(handles.axes3);
     cla;
     plot(D1);
     axes(handles.axes2);
     cla;
     plot(M);
 end
 
    case 3
        [fname,fpath]=uigetfile('*.wav');    %selection du fichier audio a ecouter
[x, fe, nbits]=wavread([fpath,fname]);
%Calculer les coefficients de prediction, signal estimée, erreur de prédiction et séquence d’autocorrélation de l’erreur de prédiction.
s= eval(get(handles.edit1,'string'));
a = lpc(x,s);
%plot(a);
est_x = filter([0 -a(2:end)],1,x);
e = x-est_x;
[acs,lags] = xcorr(e,'coeff');
%representation des resultats
% plot(x,'y');
% hold on 
% plot(est_x,'r');
plot(1:97,x(4001:4097),1:97,est_x(4001:4097),'r--'), grid
title 'signal original vs. estimation LPC'
xlabel 'nombre d'echantillons', 'ylabel 'Amplitude'
legend('signal original','signal estimé par LPC')
% %Tracer l’autocorrélation de l’erreur de prédiction.
% figure;
% plot(lags,acs), grid
% title 'Autocorrelation of the Prediction Error'
% xlabel 'Lags', ylabel 'Normalized value'
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

set(hObject, 'String', {'compression DCT', 'graphe DTW', 'représentation LPC'});

% --- Executes on selection change in popupmenu3.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes on button press in aide.
function aide_Callback(hObject, eventdata, handles)
% hObject    handle to aide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dos('notepad  fichier.txt &');



% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)

% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on key press with focus on radiobutton2 and none of its controls.
function radiobutton2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over tracer.
function tracer_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to tracer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
