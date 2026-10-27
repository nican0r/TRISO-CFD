function commandwindow(varargin)
%COMMANDWINDOW  No-op shim for headless (matlab -batch) runs.
%
% PIVlab's piv.piv_FFTmulti calls the built-in COMMANDWINDOW to bring the
% desktop command window to the front when its interpass validation discards
% most of the vectors of a pass (which the C26 spout view does by design:
% ~70-80 % of its windows hold no particles).  The built-in errors without a
% desktop ("No Desktop found") and aborts the pair; the PIV result itself is
% not affected (the warning branch only prints and beeps).  piv_c26_pivlab_run.m
% adds this folder to the path only when no desktop is running
% (~usejava('desktop')), so the GUI is unaffected.
end
