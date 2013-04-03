function varargout = fitprfgui(varargin)

% function fitprfgui
%
% this function implements a GUI for fitprf.m.  the purpose is:
% (1) to make calling the fitprf.m command a little easier.
%     this is achieved, in part, by suggesting some reasonable values
%     for the various input parameters of fitprf.m.
% (2) to include some domain-specific knowledge and defaults that
%     fitprf.m does not already have
%
% note that fitprfgui.m does not afford any functionality that
% is not already present in fitprf.m.  (in fact, the functionality
% of fitprf.m is a superset of that of fitprfgui.m.)
%
% also, note that we attempt to make the functionality of fitprfgui.m
% as transparent as possible to the user.  this is in order to avoid
% the problem of GUIs becoming impenetrable black boxes.
%
% calling fitprfgui.m will bring up a GUI (multiple instances of the GUI
% are okay).  the top part of the GUI (up through "accuracy metric:") 
% contains various GUI elements that can be changed by the user.
% these GUI elements control the various input parameters that are passed
% to fitprf.m:
%
% 1. 'stimulus': enter the name of a variable that conforms to the
%                format of the <stimulus> input of fitprf.m
% 2. 'response': enter the name of a variable that conforms to the
%                format of the <response> input of fitprf.m
% 3. 'TR': enter the time in seconds that separates successive data points
% 4. 'PRF model':
%      'GLM' means to interpret the stimulus as a standard-GLM coding of the
%        experimental design (e.g. ones indicate trial onsets, with zeros
%        elsewhere).  if this option is selected, we use a PRF model that
%        consists of delta basis functions (specifically, the '0' case of 
%        the <prfmodel> input in fitprf.m), which basically means that the PRF
%        model doesn't do anything.
%      '2D Gaussian' means to interpret the stimulus as a coding of a retinotopic
%        mapping design.  the stimulus should be time x pixels*pixels with values 
%        indicating where the aperture is located at each time point (e.g. values 
%        could be 0 or 1 with 1s indicating the location of the aperture).  if this 
%        option is selected, we use a PRF model that is an isotropic 2D Gaussian 
%        function parameterized by x (column index of peak), y (row index of peak),
%        and sigma (std dev of Gaussian); see constructprfmodel2Dgaussian.m for
%        details.
% 5. 'HRF model':
%      'delta basis' means to use delta functions to model the HRF.
%                    see constructhrfmodeldelta.m for details.
%      'DCT basis' means to use DCT-II functions to model the HRF.  by default, we
%                  include only functions that have a frequency that is less than or
%                  equal to 0.2 cycles per second.  see constructhrfmodeldct.m
%                  for details.
%      'spline' means to use a spline-based model of the HRF.  this model incorporates
%               certain default assumptions; see constructhrfmodelspline.m for details.
% 6. 'HRF duration': enter the desired duration of the HRF.  this field is editable
%      only if the 'HRF model' is the delta or DCT case.  in the spline case, the HRF
%      duration is fixed to 50 seconds.
% 7. 'separable': check this box if you want to enforce separability between the PRF 
%      and the HRF, i.e. enforce the constraint that each PRF component has the same
%      HRF up to a scale factor.  in most cases, you probably want to enforce separability.
%      cases where you would not want to enforce separability include cases where you expect 
%      differences in the HRF across trial types (e.g. due to different trial durations) and 
%      cases where you have only one trial type, in which case separability and inseparability 
%      are equivalent and we might as well not enforce separability to get faster execution.
%      note that the separability option matters only in the case where both the PRF and HRF
%      are modeled using basis functions; in the other cases, separability is implied.
% 8. 'max poly degree': select the desired maximum degree of polynomial to use for modeling
%      low-frequency drift in the time-series.
% 9. 'nonlinearity':
%      'none' means do not include any nonlinear operation.  note that when the PRF model
%             is the 'GLM' case, this is the only option that can be selected.
% 10. 'max iterations': enter the maximum number of iterations to allow in the cases
%       where we use nonlinear optimization to fit the model.
% 11. 'resampling':
%       'none' means just fit the model to all the data
%       'n-fold cross-validation' means systematically cross-validate on each run.
%         on the first fit, we use all runs except the first to fit and we then cross-validate
%         on the first run.  on the second fit, we use all runs except the second to fit and
%         we then cross-validate on the second run.  and so forth.  the point of cross-validation
%         is to obtain an unbiased measure of model accuracy.
%       'bootstrap (30 reps)' means to perform 30 separate model fits.  for each model fit,
%         we draw, with replacement, from the available set of runs and fit the model parameters
%         to these data.  the point of bootstrapping is obtain estimates of the standard error
%         on the various model parameters.  (note that we bootstrap runs, not data points.)
% 12. 'tolerance': enter the tolerance to use in the cases where we use nonlinear
%       optimization to fit the model.
% 13. 'large-scale': check this box if you want to use the 'LargeScale' option in the cases
%       where we use nonlinear optimization to fit the model.  whether or not you check this box
%       will affect the speed of convergence and susceptibility to local minima.  based on
%       my observations, the optimization may be much faster if you leave the box unchecked.
% 14. 'observe fitting': check this box if you want to observe the progress of the nonlinear
%       optimization in a figure window.  see outputfcnplot.m for more details (we plot 
%       after every iteration).
% 15. 'HRF normalization': when fitting a model that allows flexibility in the shape of the HRF,
%                          it is not immediately obvious how to associate a single beta value
%                          that represents the amplitude of the response to a given PRF component.
%                          this is because scaling the HRF by a constant has the same effect as
%                          dividing PRF components by that same constant.
%       'standard' means to use a normalization strategy wherein we force the estimated HRF to
%         have a peak of exactly one.  we do this by estimating the peak (which can be positive
%         or negative) of the HRF using calchrfpeak.m and then dividing the HRF by this estimated
%         peak.  note that we use certain defaults in the call to calchrfpeak.m.
% 16. 'derive beta weights': check this box if, after the fitting of the model, you would like 
%       to perform some convenient post-hoc calculations involving the derivation of beta weights.
%       it is recommended that you leave this box checked since the derivation of beta weights does
%       not take much computational time.  when the resampling type is n-fold cross-validation,
%       there is an additional pop-up menu:
%         'group' means to derive beta weights for training runs as a group and testing runs as a group
%                 (thus, there is a single beta weight for training and a single weight for testing)
%         'individual' means to derive beta weights for individual runs (thus, there may be multiple
%                      beta weights for training and multiple beta weights for testing)
%         'both' means to do 'group' and also do 'individual'
%       for additional details on the group and individual cases, see fitprf.m.
% 17. 'accuracy metric': this controls how model accuracy is quantified.  note that polynomials are 
%       projected out from both the data and the model before calculation of model accuracy.  also, 
%       note that in the cross-validation case, we quantify model accuracy on runs that are not used
%       to the fit model.  see fitprf.m for more details on exactly how model accuracy is calculated.
%         'R^2' means to quantify model accuracy as percent variance of the data explained by the model.
%
% the bottom part of the GUI contain buttons that actually cause some action to occur:
%
% 1. 'Run and Save to': this button causes the fitprf.m command to be issued, with the outputs
%      encapsulated inside a struct.  this struct is saved to the workspace to a variable
%      whose name is determined by the adjoining text box.  the command that does all of this
%      is reported to the command window.
% 2. 'Run and Save to Workspace': this causes the fitprf.m command to be issued, with the outputs
%      being saved directly to the workspace.  (the names of the variables used are those 
%      described in the documentation of fitprf.m.)  the command that does all of this
%      is reported to the command window.
% 3. 'Get Command': this merely reports the command that would be run if you were to click the
%      'Run and Save to Workspace' button.
%
% history:
% 2010/12/06 - implement 2D Gaussian PRF model and clean up GUI dependencies
% 2010/12/06 - initial version
%
% some technical details:
% - certain inputs of fitprf.m are not explored by fitprfgui.m.  these include:
%   - <extraregressors> is always []
%   - <flag> is never {2 X Y}
%   - <ar> is always []

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FITPRFGUI M-file for fitprfgui.fig
%      FITPRFGUI, by itself, creates a new FITPRFGUI or raises the existing
%      singleton*.
%
%      H = FITPRFGUI returns the handle to a new FITPRFGUI or the handle to
%      the existing singleton*.
%
%      FITPRFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITPRFGUI.M with the given input arguments.
%
%      FITPRFGUI('Property','Value',...) creates a new FITPRFGUI or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fitprfgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fitprfgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fitprfgui

% Last Modified by GUIDE v2.5 05-Dec-2010 11:27:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fitprfgui_OpeningFcn, ...
                   'gui_OutputFcn',  @fitprfgui_OutputFcn, ...
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

% --- Executes just before fitprfgui is made visible.
function fitprfgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fitprfgui (see VARARGIN)

% Choose default command line output for fitprfgui
handles.output = hObject;

%%%%%

% we need to store the values of some fields since we do some basic error-checking
handles.trval = 1;
handles.hrfdurationval = 50;
handles.maxiterationsval = Inf;
handles.toleranceval = 1e-6;
  % this isn't really necessary since we edited the .fig file to have these
  % values already, but let's do it anyway.
set(handles.tr,'String',num2str(handles.trval));
set(handles.hrfduration,'String',num2str(handles.hrfdurationval));
set(handles.maxiterations,'String',num2str(handles.maxiterationsval));
set(handles.tolerance,'String',num2str(handles.toleranceval));

% some of the GUI elements imply certain dependencies depending on their
% values.  so, let's take care of those dependencies now.
handledependencies(handles);

%%%%%

% Update handles structure
guidata(handles.output, handles);

% UIWAIT makes fitprfgui wait for user response (see UIRESUME)
% uiwait(handles.fitprfguifig);


% --- Outputs from this function are returned to the command line.
function varargout = fitprfgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function tr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tr_Callback(hObject, eventdata, handles)
% hObject    handle to tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(hObject,'String'));
if ~isfinite(val) || val <= 0
  set(hObject,'String',num2str(handles.trval));
else
  handles.trval = val;
end
guidata(handles.output,handles);



% --- Executes on selection change in prfmodel.
function prfmodel_Callback(hObject, eventdata, handles)
% hObject    handle to prfmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns prfmodel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from prfmodel

handledependencies(handles);



% --- Executes during object creation, after setting all properties.
function prfmodel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prfmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in hrfmodel.
function hrfmodel_Callback(hObject, eventdata, handles)
% hObject    handle to hrfmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns hrfmodel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from hrfmodel

handledependencies(handles);



% --- Executes during object creation, after setting all properties.
function hrfmodel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hrfmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hrfduration_Callback(hObject, eventdata, handles)
% hObject    handle to hrfduration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hrfduration as text
%        str2double(get(hObject,'String')) returns contents of hrfduration as a double

val = str2double(get(hObject,'String'));
if ~isfinite(val) || val <= 0
  set(hObject,'String',num2str(handles.hrfdurationval));
else
  handles.hrfdurationval = val;
end
guidata(handles.output,handles);


% --- Executes during object creation, after setting all properties.
function hrfduration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hrfduration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in separable.
function separable_Callback(hObject, eventdata, handles)
% hObject    handle to separable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of separable

handledependencies(handles);



% --- Executes on selection change in nonlinearity.
function nonlinearity_Callback(hObject, eventdata, handles)
% hObject    handle to nonlinearity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nonlinearity contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nonlinearity


% --- Executes during object creation, after setting all properties.
function nonlinearity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonlinearity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxiterations_Callback(hObject, eventdata, handles)
% hObject    handle to maxiterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxiterations as text
%        str2double(get(hObject,'String')) returns contents of maxiterations as a double

val = str2double(get(hObject,'String'));
if isnan(val) || val < 1
  set(hObject,'String',num2str(handles.maxiterationsval));
else
  handles.maxiterationsval = round(val);
end
guidata(handles.output,handles);


% --- Executes during object creation, after setting all properties.
function maxiterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxiterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in resampling.
function resampling_Callback(hObject, eventdata, handles)
% hObject    handle to resampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns resampling contents as cell array
%        contents{get(hObject,'Value')} returns selected item from resampling

handledependencies(handles);



% --- Executes during object creation, after setting all properties.
function resampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tolerance as text
%        str2double(get(hObject,'String')) returns contents of tolerance as a double

val = str2double(get(hObject,'String'));
if ~isfinite(val) || val <= 0
  set(hObject,'String',num2str(handles.toleranceval));
else
  handles.toleranceval = val;
end
guidata(handles.output,handles);


% --- Executes during object creation, after setting all properties.
function tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in largescale.
function largescale_Callback(hObject, eventdata, handles)
% hObject    handle to largescale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of largescale


% --- Executes on button press in observefitting.
function observefitting_Callback(hObject, eventdata, handles)
% hObject    handle to observefitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of observefitting


% --- Executes on selection change in hrfnormalization.
function hrfnormalization_Callback(hObject, eventdata, handles)
% hObject    handle to hrfnormalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns hrfnormalization contents as cell array
%        contents{get(hObject,'Value')} returns selected item from hrfnormalization


% --- Executes during object creation, after setting all properties.
function hrfnormalization_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hrfnormalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function stimulus_Callback(hObject, eventdata, handles)
% hObject    handle to stimulus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stimulus as text
%        str2double(get(hObject,'String')) returns contents of stimulus as a double


% --- Executes during object creation, after setting all properties.
function stimulus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimulus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function response_Callback(hObject, eventdata, handles)
% hObject    handle to response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of response as text
%        str2double(get(hObject,'String')) returns contents of response as a double


% --- Executes during object creation, after setting all properties.
function response_CreateFcn(hObject, eventdata, handles)
% hObject    handle to response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function runandsavevar_Callback(hObject, eventdata, handles)
% hObject    handle to runandsavevar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runandsavevar as text
%        str2double(get(hObject,'String')) returns contents of runandsavevar as a double



% --- Executes during object creation, after setting all properties.
function runandsavevar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runandsavevar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in accuracymetric.
function accuracymetric_Callback(hObject, eventdata, handles)
% hObject    handle to accuracymetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns accuracymetric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from accuracymetric


% --- Executes during object creation, after setting all properties.
function accuracymetric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to accuracymetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in derivebetaweights.
function derivebetaweights_Callback(hObject, eventdata, handles)
% hObject    handle to derivebetaweights (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of derivebetaweights

handledependencies(handles);



% --- Executes on selection change in derivebetaweightstype.
function derivebetaweightstype_Callback(hObject, eventdata, handles)
% hObject    handle to derivebetaweightstype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns derivebetaweightstype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from derivebetaweightstype


% --- Executes during object creation, after setting all properties.
function derivebetaweightstype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to derivebetaweightstype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in maxpolydeg.
function maxpolydeg_Callback(hObject, eventdata, handles)
% hObject    handle to maxpolydeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns maxpolydeg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from maxpolydeg


% --- Executes during object creation, after setting all properties.
function maxpolydeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxpolydeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function separable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to separable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in runandsave.
function runandsave_Callback(hObject, eventdata, handles)
% hObject    handle to runandsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

runit(get(handles.runandsavevar,'String'),1,handles);

% --- Executes on button press in runandsavetoworkspace.
function runandsavetoworkspace_Callback(hObject, eventdata, handles)
% hObject    handle to runandsavetoworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

runit('[params,paramsse,r,rrun,polyparams,polymeans,numiters,hrf,betas,signal,drift]',1,handles);

% --- Executes on button press in getcommand.
function getcommand_Callback(hObject, eventdata, handles)
% hObject    handle to getcommand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

runit('[params,paramsse,r,rrun,polyparams,polymeans,numiters,hrf,betas,signal,drift]',0,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handledependencies(handles)

% prfmodel and resampling and derivebetaweights affect derivebetaweights and derivebetaweightstype
if get(handles.prfmodel,'Value')==2  % if 2D Gaussian case, turn everything off
  set(handles.derivebetaweightstext,'Enable','off');
  set(handles.derivebetaweights,'Value',0);
  set(handles.derivebetaweights,'Enable','off');
  set(handles.derivebetaweightstype,'Enable','off');
else
  set(handles.derivebetaweightstext,'Enable','on');
  set(handles.derivebetaweights,'Enable','on');
  switch get(handles.resampling,'Value')
  case {1 3}  % if no resampling or bootstrapping, then the type doesn't apply
    set(handles.derivebetaweightstype,'Value',1);
    set(handles.derivebetaweightstype,'Enable','off');
  case 2  % if cross-validation, then the type does apply (if the checkbox is checked)
    if get(handles.derivebetaweights,'Value')
      set(handles.derivebetaweightstype,'Enable','on');
    else
      set(handles.derivebetaweightstype,'Enable','off');
    end
  end
end

% prfmodel affects nonlinearity
switch get(handles.prfmodel,'Value')
case {1 2}  % if GLM or 2D-Gaussian model, no nonlinearity!
  set(handles.nonlinearity,'Value',1);
  set(handles.nonlinearity,'Enable','off');
end

% hrfmodel affects hrfduration
switch get(handles.hrfmodel,'Value')
case {1 2}  % if delta or DCT, we want the duration
  set(handles.hrfdurationtext,'Enable','on');
  set(handles.hrfduration,'Enable','on');
  set(handles.hrfdurationtext2,'Enable','on');
case 3  % if spline, we don't want the duration [and we act like it's 50]
  handles.hrfdurationval = 50;
  set(handles.hrfduration,'String',num2str(handles.hrfdurationval));
  set(handles.hrfdurationtext,'Enable','off');
  set(handles.hrfduration,'Enable','off');
  set(handles.hrfdurationtext2,'Enable','off');
end

% prfmodel and hrfmodel affect separable
if get(handles.prfmodel,'Value')==1 && ismember(get(handles.hrfmodel,'Value'),[1 2])
  set(handles.separabletext,'Enable','on');
  set(handles.separable,'Enable','on');
else
  set(handles.separable,'Value',1);
  set(handles.separabletext,'Enable','off');
  set(handles.separable,'Enable','off');
end

% save
guidata(handles.output,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runit(str,dorun,handles)

% <str> is the output argument to assign to
% <dorun> is whether to actually run the command

% PRFMODEL
switch get(handles.prfmodel,'Value')
case 1
  prfmodel = '0';
case 2
  prfmodel = sprintf('constructprfmodel2Dgaussian(cellfunfirst(@(x)sqrt(size(x,2)),%s))',get(handles.stimulus,'String'));
end

% HRFMODEL
switch get(handles.hrfmodel,'Value')
case 1  % delta
  hrfmodel = sprintf('constructhrfmodeldelta(%0.10g,%0.10g)',handles.hrfdurationval,handles.trval);
case 2  % dct
  hrfmodel = sprintf('constructhrfmodeldct(%0.10g,%0.10g)',handles.hrfdurationval,handles.trval);
case 3  % spline
  hrfmodel = sprintf('constructhrfmodelspline(50,%0.10g,[2.5:2.5:20 30 40],[.5 1 1 .5 zeros(1,6)])',handles.trval);
end

% MODE
switch get(handles.nonlinearity,'Value')
case 1
  mode = '0';
end

% WANTRESAMPLE
switch get(handles.resampling,'Value')
case 1
  wantresample = '0';
case 2
  wantresample = '''n-fold''';
case 3
  wantresample = '{-30 [] []}';
end

% EXTRAOPT
if get(handles.largescale,'Value')
  extraopt = '{''LargeScale'' ''on''';
else
  extraopt = '{''LargeScale'' ''off''';
end
if get(handles.observefitting,'Value')
  extraopt = [extraopt ' ''OutputFcn'' @(a,b,c)outputfcnplot(a,b,c,1)}'];
else
  extraopt = [extraopt '}'];
end

% HRFNORMFUN
switch get(handles.hrfnormalization,'Value')
case 1
  hrfnormfun = sprintf('@(x)calchrfpeak(x,%0.10g)',handles.trval);
end

% DERIVEMODE
if get(handles.derivebetaweights,'Value')
  if ismember(get(handles.resampling,'Value'),[1 3])
    derivemode = '1';
  else  % ok, it must be 2
    switch get(handles.derivebetaweightstype,'Value')
    case 1  % group
      derivemode = '[1 1]';
    case 2  % individual
      derivemode = '[2 2]';
    case 3  % both
      derivemode = '[3 3]';
    end
  end
else
  derivemode = '[]';
end

% METRIC
switch get(handles.accuracymetric,'Value')
case 1
  metric = '[]';
end

% do it
expr = sprintf('%s = fitprf(%s,%s,%s,%s,%d,%d,[],%s,%d,%s,%0.10g,%s,[],%s,%s,%s);', ...
  str,get(handles.stimulus,'String'),get(handles.response,'String'), ...
  prfmodel,hrfmodel,get(handles.separable,'Value'),get(handles.maxpolydeg,'Value')-1,mode, ...
  handles.maxiterationsval,wantresample,handles.toleranceval,extraopt,hrfnormfun,derivemode,metric);
if dorun
  fprintf(['we are running the following command:\n' expr '\n']);
  evalin('base',expr);
else
  fprintf(['the command is:\n' expr '\n']);
end
