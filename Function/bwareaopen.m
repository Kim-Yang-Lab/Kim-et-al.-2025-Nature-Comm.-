function bw2 = bwareaopen(varargin)
%BWAREAOPEN Remove small objects from binary image.
%   BW2 = BWAREAOPEN(BW,P) removes from a binary image all connected
%   components (objects) that have fewer than P pixels, producing another
%   binary image BW2.  This operation is known as an area opening.  The
%   default connectivity is 8 for two dimensions, 26 for three dimensions,
%   and CONNDEF(NDIMS(BW),'maximal') for higher dimensions. 
%
%   BW2 = BWAREAOPEN(BW,P,CONN) specifies the desired connectivity.  CONN
%   may have the following scalar values:  
%
%       4     two-dimensional four-connected neighborhood
%       8     two-dimensional eight-connected neighborhood
%       6     three-dimensional six-connected neighborhood
%       18    three-dimensional 18-connected neighborhood
%       26    three-dimensional 26-connected neighborhood
%
%   Connectivity may be defined in a more general way for any dimension by
%   using for CONN a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The
%   1-valued elements define neighborhood locations relative to the center
%   element of CONN.  CONN must be symmetric about its center element.
%
%   Class Support
%   -------------
%   BW can be a logical or numeric array of any dimension, and it must be
%   nonsparse.
%
%   BW2 is logical.
%
%   Example
%   -------
%   % Remove all objects in the image text.png containing fewer than 50
%   % pixels.
%
%       BW = imread('text.png');
%       BW2 = bwareaopen(BW,50);
%       imshow(BW);
%       figure
%       imshow(BW2)
%
%   See also BWCONNCOMP, CONNDEF, REGIONPROPS.

%   Copyright 1993-2019 The MathWorks, Inc.

% Input/output specs
% ------------------
% BW:    N-D real full matrix
%        any numeric class
%        sparse not allowed
%        anything that's not logical is converted first using
%          bw = BW ~= 0
%        Empty ok
%        Inf's ok, treated as 1
%        NaN's ok, treated as 1
%
% P:     double scalar
%        nonnegative integer
%
% CONN:  connectivity
%
% BW2:   logical, same size as BW
%        contains only 0s and 1s.

matlab.images.internal.errorIfgpuArray(varargin{:});
[bw,p,conn] = parse_inputs(varargin{:});

CC = bwconncomp(bw,conn);
area = cellfun(@numel, CC.PixelIdxList);

idxToKeep = CC.PixelIdxList(area >= p);
idxToKeep = vertcat(idxToKeep{:});

bw2 = false(size(bw));
bw2(idxToKeep) = true;


%%%
%%% parse_inputs
%%%
function [bw,p,conn] = parse_inputs(varargin)

narginchk(2,3)

bw = varargin{1};
validateattributes(bw,{'numeric', 'logical'},{'real', 'nonsparse'},mfilename,'BW',1);
if ~islogical(bw)
    bw = bw ~= 0;
end

p = varargin{2};
validateattributes(p,{'double'},{'scalar' 'integer' 'nonnegative'},...
    mfilename,'P',2);

if (nargin >= 3)
    conn = varargin{3};
else
    conn = conndef(ndims(bw),'maximal');
end
iptcheckconn(conn,mfilename,'CONN',3)
