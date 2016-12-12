function [delH_NN, delS_NN, delG_NN, ...
	  delH_AT_closing_penalty, delS_AT_closing_penalty, ...
	  delG_AT_closing_penalty,...
	  delH_mismatch, delS_mismatch, delG_mismatch, ...
	  delH_init, delS_init ...
	 ] = get_NN_parameters( T );

if ~exist( 'T' )
  T = 273.15 + 37;
end

%delH = -1 * [ 9.1, 6.5, 7.8, 8.6; ...
%	      5.8, 11.0, 11.9, 7.8; ...
%	      5.6, 11.1, 11.0, 6.5; ...
%	      6.0, 5.6, 5.8, 9.1];
%delS = -1 * [ 24.0, 17.3, 20.8, 23.9; ...
%	      12.9, 26.6, 27.8, 20.8; ...
%	      13.5, 26.7, 26.6, 17.3; ...
%	      16.9, 13.5, 12.9, 24.0];

% From SantaLucia, Jr, PNAS, 1998.
%delH = [ -7.9, -8.4,  -7.8, -7.2; ...
%	 -8.5, -8.0, -10.6, -7.8; ...
%	 -8.2, -9.8,  -8.0, -8.4; ...
%	 -7.2, -8.2,  -8.5, -7.9];
%delS = [ -22.2, -22.4, -21.0, -20.4; ...
%	 -22.7, -19.9, -27.2, -21.0; ...
%	 -22.2, -24.4, -19.9, -22.4; ...
%	 -21.3, -22.2, -22.7, -22.2];

% From SantaLucia, Jr, Ann Rev 2004
% A C G T
delH_NN = [ -7.6, -8.4,  -7.8, -7.2; ...
	 -8.5, -8.0,  -9.8, -7.8; ...
	 -8.2, -9.8,  -8.0, -8.4; ...
	 -7.2, -8.2,  -8.5, -7.6];
delS_NN = [ -21.3, -22.4, -21.0, -20.4; ...
	    -22.7, -19.9, -24.4, -21.0; ...
	    -22.2, -24.4, -19.9, -22.4; ...
	    -21.3, -22.2, -22.7, -21.3];

delG_NN = delH_NN - (T * delS_NN)/1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From SantaLucia, Jr, Ann Rev 2004
delH_AT_closing_penalty = [ 2.2 0 0 2.2];
delS_AT_closing_penalty = [ 6.9 0 0 6.9];
delG_AT_closing_penalty = delH_AT_closing_penalty - (T * delS_AT_closing_penalty)/1000;



% Following also from Santaucia/Hicks.
% They didn't compile delH and delS, so I should probably go
% back to the original references. Pain in the ass!

% AX/TY
delH_mismatch(:,:,1) = [...
     1.2 2.3 -0.6 -7.6;
     5.3 0.0 -8.4 0.7; ...
    -0.7 -7.8 -3.1 1.0; ...
    -7.2 -1.2 -2.5 -2.7 ] ;

%  CX/GY
delH_mismatch(:,:,2) = [...
    -0.9 1.9 -0.7 -8.5;...
    0.6 -1.5 -8.0 -0.8; ...
    -4.0 -10.6 -4.9 -4.1; ...
    -7.8 -1.5 -2.8 -5.0 ] ;
    
%  GX/CY
delH_mismatch(:,:,3) = [...
    -2.9 5.2 -0.6 -8.2; ...
    -0.7 3.6 -9.8 2.3; ...
    0.5 -8.0 -6.0 3.3; ...
    -8.4 5.2 -4.4 -2.2 ];
    
%  TX/AY
delH_mismatch(:,:,4) = [...
    4.7 3.4 0.7 -7.2; ...
    7.6 6.1 -8.2 1.2; ...
    3.0 -8.5 1.6 -0.1; ...
    -7.6 1.0 -1.3 0.2 ];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AX/TY
delS_mismatch(:,:,1) = [...
     1.7  4.6 -2.3 -21.3; ...
     14.6 -4.4 -22.4 0.2; ...
    -2.3 -21.0 -9.5 0.9; ...
    -20.4 -6.2 -8.3 -10.8 ];
    
%  CX/GY
delS_mismatch(:,:,2) = [...
    -4.2 3.7 -2.3 -22.7; ...
    -0.6 -7.2 -19.9 -4.5; ...
    -13.2 -27.2 -15.3 -11.7; ...
    -21.0 -6.1 -8.0 -15.8];
    
%  GX/CY
delS_mismatch(:,:,3) = [...
    -9.8 14.2 -1.0 -22.2; ...
    -3.8 8.9 -24.4 5.4; ...
    3.2 -19.9 -15.8 10.4; ...
    -22.4 13.5 -12.3 -8.4 ];
    
%  TX/AY
delS_mismatch(:,:,4) = [...
    12.9 8.0 0.7 -21.3; ...
    20.2 16.4 -22.2 0.7; ...
    7.4 -22.7 3.6 -1.7; ...
    -21.3 0.7 -5.3 -1.5 ] ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AX/TY
delG_mismatch(:,:,1) = [...
    0.61   0.88  0.14 -1.0; ...
    0.77   1.33  -1.44  0.64; ...
    0.02  -1.28 -0.13  0.71; ...
    -0.88 0.73   0.07  0.69 ];

%  CX/GY
delG_mismatch(:,:,2) = [...
    0.43 0.75 0.03 -1.45; ...
    0.79 0.70 -1.84 0.62; ...
    0.11 -2.17 -0.11 -0.47; ...
    -1.28 0.40 -0.32 -0.13 ];

%  GX/CY
delG_mismatch(:,:,3) = [...
    0.17 0.81 -0.25 -1.30;...
    0.47 0.79 -2.24 0.62; ...
    -0.52 -1.84 -1.11 0.08; ...
    -1.44 0.98 -0.59 0.45];

%  TX/AY
delG_mismatch(:,:,4) = [...
    0.69 0.92 0.42 -0.58;...
    1.33 1.05 -1.30 0.97; ...
    0.74 -1.45 0.44 0.43; ...
    -1.00 0.75 0.34 0.68 ];

delH_init = 0.2;
delS_init = -5.7;


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For debugging...
delG_mismatch_37 = delH_mismatch - (273.15+37)*delS_mismatch/1000;

blah1 =  reshape( delG_mismatch_37, 1, 64);
blah2 =  reshape( delG_mismatch ,  1, 64);
pos = find( abs( blah1 - blah2)> 0.2 & blah2 < 999);
blah2( pos )
blah1(pos) - blah2(pos)

plot(blah1, blah2,'.' );
axis([-20 20 -20 20]);

