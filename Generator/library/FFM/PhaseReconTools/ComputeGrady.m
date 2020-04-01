function phase_grady = ComputeGrady(phase, pupil_mask)
%
%     phase_grady = ComputeGrady(phase, pupil_mask);
%
%  This function can be used to compute gradients of a phase, assuming
%  a Fried geometry of the WFS.
%
%  Input:
%    phase      - array containing phases -- could be 3D
%
%  Optional Input:
%    pupil_mask - array defining pupil mask -- could be 3D.
%                     2D array ==> pupil mask is the same for each phase
%                     3D array ==> size(pupil_mask) = size(phase)
%                 If this is given, then the phase is padded outside
%                 the pupil mask using an inpainting technique.
%
%  Output:
%    phase_gradx - computed gradients
%                  array same size as phase
%
%
if nargin == 1
  pupil_mask = [];
end

[m_ap, n_ap, n_frames] = size(phase);
if m_ap ~= n_ap
  error('Expected input frames to be square')
end

Dy = GradyMatrix(n_ap);

if ~isempty(pupil_mask)
  if size(pupil_mask,3) ~= size(phase,3)
    pupil_mask = pupil_mask(:,:,ones(1,n_frames));
  end
  for k = 1:n_frames
    phase(:,:,k) = MosaicPadding(phase(:,:,k), pupil_mask(:,:,k));
  end
end

phase_grady = zeros(size(phase));
for k = 1:n_frames
  phase_grady(:,:,k) = pupil_mask(:,:,k).*reshape(Dy*reshape(phase(:,:,k),n_ap*n_ap,1),n_ap,n_ap);
end

