addpath(genpath('lib'))

main03_ir('ir_type', 'square');
close all
main03_ir('ir_type', 'erp');
close all

main04_snr;
close all

main05_emph('ir_type', 'square');
close all
main05_emph('ir_type', 'erp2');
close all

main06_scale;
close all

main07_shift;
close all

main08_coch;
close all

main09_eeg;
close all

main10_jitter;
close all