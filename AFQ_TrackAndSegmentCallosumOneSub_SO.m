function AFQ_TrackAndSegmentCallosumOneSub_SO(dtDir, varargin)
% Track, segment and clean fibers for a single subject
%
% AFQ_TrackAndSegment(dtDir, [params])
%
% This function will track, segment and clean fibers. The resulting files
% will be save so that future runs of AFQ will be faster
%
% Inputs:
%
% dtDir    = Path to the directory containing the dt6 file
% params   = Parameters to change in the afq structure. See AFQ_Create
%
% Copyright Jason D. Yeatman, December 2012

sub_dirs{1} = dtDir;
afq = AFQ_Create('sub_dirs',sub_dirs,'sub_group',1);
% Set to overwrite previous fibers
afq=AFQ_set(afq,'overwritefibers');
% Set coupute norms to false
afq.params.computenorms =0;

% add fake values in to convince AFQ_run not to compute anything
%afq.vals.fa{1}=zeros(1,100);
if nargin>1
    for ii = 1:2:length(varargin)
        afq = AFQ_set(afq,varargin{ii},varargin{ii+1});
    end
end

% run AFQ_SegmentCallosum
AFQ_SegmentCallosum(afq);