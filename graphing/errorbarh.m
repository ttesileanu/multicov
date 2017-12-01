function hh = errorbarh(varargin)
% ERRORBARH Draw horizontal error bars.
%   ERRORBARH(X, Y, L, R) plots horizontal error bars for Y vs. X with
%   errors given by L and R. More precisely, the error bar for point i goes
%   from X(i) - L(i) to X(i) + R(i). The inputs can be vectors or matrices;
%   when they are matrices, each column produces a separate set of error
%   bars.
%
%   Note that unlike Matlab's errorbar function, this does not plot the
%   data itself; it only draws the error bars. Because of this, if a plot
%   is open, this will add the error bars to that plot instead of opening a
%   new one (even if hold wasn't on).
%
%   ERRORBARH(X, Y, E) uses symmetric errorbars (L = R = E).
%
%   ERRORBARH(Y, E, ...) uses 1:size(Y, 1) as x coordinates.
%
%   Additional parameters passed to ERRORBARH can be used to modify the
%   color and style of the error bars. See Matlab's plot documentation for
%   details.
%
%   H = ERRORBARH(...) returns a vector of line handles.
%
%   Example:
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      errorbarh(x, y, e);
%   draws symmetric horizontal error bars of unit standard deviation.
%
% This code is modified from Jos van der Geest's Mathworks submission
% http://www.mathworks.com/matlabcentral/fileexchange/3963-herrorbar.
%
% See also: ERRORBAR, PLOT.

% Tiberiu Tesileanu (2014)

%   Jos van der Geest
%   email: jos@jasen.nl
%
%   File history:
%   August 2006 (Jos): I have taken back ownership. I like to thank Greg Aloe from
%   The MathWorks who originally introduced this piece of code to the
%   Matlab File Exchange. 
%   September 2003 (Greg Aloe): This code was originally provided by Jos
%   from the newsgroup comp.soft-sys.matlab:
%   http://newsreader.mathworks.com/WebX?50@118.fdnxaJz9btF^1@.eea3ff9
%   After unsuccessfully attempting to contact the orignal author, I
%   decided to take ownership so that others could benefit from finding it
%   on the MATLAB Central File Exchange.

% Copyright (c) 2009, Jos van der Geest
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% line-style arguments start with the first string argument
i = find(cellfun(@ischar, varargin), 1);
linestyle_args = varargin(i:end);
if ~isempty(i)
    numeric_args = varargin(1:i-1);
else
    numeric_args = varargin;
end
if length(numeric_args) == 2        % y, e
    y = numeric_args{1};
    l = numeric_args{2};
    r = l;
    % need to invent x coordinates
    x = repmat((1:size(y, 1))', 1, size(y, 2));
elseif length(numeric_args) > 2
    x = numeric_args{1};
    y = numeric_args{2};
    l = numeric_args{3};
    if length(numeric_args) == 3    % x, y, e
        r = l;
    elseif length(numeric_args) == 4
        r = numeric_args{4};
    else
        error([mfilename ':toomany'], 'Too many arguments.');
    end
else
    error([mfilename ':toofew'], 'Need at least two arguments.');
end

% make sure errors are positive
l = abs(l);
r = abs(r);

% make sure they are numeric (but logical is fine, too)
checkn = @(v) isnumeric(v) || islogical(v);
if ~checkn(x) || ~checkn(y) || ~checkn(l) || ~checkn(r)
    error([mfilename ':nonnum'], 'Arguments must be numeric.');
end

if ~isequal(size(x), size(y)) || ~isequal(size(x), size(l)) || ~isequal(size(x), size(r))
    error([mfilename ':badsize'], 'The sizes of X, Y, L and U must be the same.');
end

% decide on the size of the tee
% this is half of it
htee = (max(y(:)) - min(y(:))) / 100;
xl = x - l;
xr = x + r;
ytop = y + htee;
ybot = y - htee;
% npt = number of points, n = number of plots
[npt, n] = size(y);

% plot the bars
hold_state = ishold;
hold on;

% build up nan-separated vector for bars
xb = zeros(npt*9, n);
xb(1:9:end, :) = xl;
xb(2:9:end, :) = xl;
xb(3:9:end, :) = NaN;
xb(4:9:end, :) = xl;
xb(5:9:end, :) = xr;
xb(6:9:end, :) = NaN;
xb(7:9:end, :) = xr;
xb(8:9:end, :) = xr;
xb(9:9:end, :) = NaN;

yb = zeros(npt*9, n);
yb(1:9:end, :) = ytop;
yb(2:9:end, :) = ybot;
yb(3:9:end, :) = NaN;
yb(4:9:end, :) = y;
yb(5:9:end, :) = y;
yb(6:9:end, :) = NaN;
yb(7:9:end, :) = ytop;
yb(8:9:end, :) = ybot;
yb(9:9:end, :) = NaN;

if isempty(linestyle_args)
    linestyle_args = {'k'};
end

h = plot(xb, yb, linestyle_args{:});

if ~hold_state
    hold off;
end

if nargout > 0
    hh = h;
end

end