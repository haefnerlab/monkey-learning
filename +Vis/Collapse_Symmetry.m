function [I, norm] = Collapse_Symmetry(I, norm, zero_idx, varargin)
%COLLAPSE_SYMMETRY with image_pref_orientation, squashes symmetries of 
% smoothed image I together. The column in I corresponding to angle 0 
% given in zero_idx
I = I .* norm;

% shift so that rotation 45 is in the center (0/90 in quadrant centers)
sz = size(I,1);
target = floor(sz / 4) + 1;
shift = zero_idx - target;

I = wshift('2D', I, shift);
norm = wshift('2D', norm, shift);

if any(strcmpi('target_swap', varargin))
    % flip across diagonal (not the transpose)
    norm = norm + fliplr(rot90(norm));
    I = I + fliplr(rot90(I));
end

if any(strcmpi('plus_minus_swap', varargin))
    % corr(x+eps,y+eps) = corr(x-eps,y-eps)
    % 1. shift to put 0-idx on the edge (90 in the middle).
    % 2. flip on non-transpose diagonal
    % 3. shift back
    shift_ = -shift + zero_idx - 1;
    I = wshift('2D', I, shift_);
    norm = wshift('2D', norm, shift_);
    
    norm = norm + fliplr(rot90(norm));
    I = I + fliplr(rot90(I));
    
    I = wshift('2D', I, -shift_);
    norm = wshift('2D', norm, -shift_);
end

% Undo shifts
I = wshift('2D', I, -shift);
norm = wshift('2D', norm, -shift);

I = I ./ norm;

end