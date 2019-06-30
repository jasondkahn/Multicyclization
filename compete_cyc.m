function [t,concs] = compete_cyc(end_time,M_input_vec,kc,kb,varargin)
% Calculate amts of products in ligation reactions given the parameters and
%    monomer input concentrations
% The optional args in varargin are kcd_vec,dfA_vec,dfB_vec
% Time should be just an end time
% M_input_vec is a vector of input monomer concentrations
% kc is a single number or else a vector of the same length as M_input_vec
%    for unimolecular  cyclization rate constants
% kb is a single bimolecular ligation rate constant for all (I suppose one could elaborate to have different
%    values for AA, AB, and BB ligations
% kcd_vec,dfA_vec,dfB_vec can be left blank, or they can be single numbers,
%    or vectors of the same length as M_input_vec
% kcd_vec is a vector of cyclization rate constants for dimer circles
% kcd_multi is the max of the kcd_vec
% dfA_vec is  a vector of the dead end fractions for the A end for each class
% similarly dfB_vec
% Returns a time vector and a matrix with the each column being the amt of
%    each monomer and product
% Example of use: 
% time = (0:120)'
%df = 0;
%kc = 0.04;
%kb = 0.0008;
%M0 = 5;
%[products] = bicycle_deadends(time,kc,kb,df,M0);
%figure
%plot(time,products,'r','LineWidth',2)
%defaults
maxnumopts = 3;
num_sets = length(M_input_vec);
% Set defaults, but recognize that the kcd is probably strongly length
% dependent in interesting length ranges
kcd_num = 3*max(kc);
dfA_num = 0;
dfB_num = dfA_num;
optargs = {kcd_num,dfA_num,dfB_num};
%
numvarargs = length(varargin);
if numvarargs > maxnumopts
    error('compete_cyc:TooManyInputs', ...
        'compete_cyc allows at most %d optional inputs',maxnumopts);
end
%
% Now put the new input values from varargin into the optargs cell array, 
% overwriting the default values.
%
% But if the user has entered '' then we do not want to replace the default
% There is probably a more elegant way to do this.
aretheythere = cellfun(@isempty,varargin);
for ivarg = 1:length(varargin)
  if ~aretheythere(ivarg)
      optargs(ivarg) = varargin(ivarg);
  end
end
%optargs(1:numvarargs) = varargin;
[kcd_num,dfA_num,dfB_num] = optargs{:};
% Allow for all inputs to be either single numbers or vectors of the
% correct length.
if (length(kc) == num_sets)
    kc_vec = kc;
elseif (length(kc) == 1)
    kc_vec = ones(1,num_sets).*kc;
else
    error('compete_cyc: kc should have length 1 or %d = number of monomer sets',num_sets);
end
% In theory there are kcd's for each pair of sets, with the number of
% interclass pairs given by te pascal triangle number below, could do
% algorithmically of course
pascalp = [0 1 3 6 10 15 21 28 36 45 55 66 78 91 105 120 136 155];
% etc. inter-pairings are triangular numbers from the Pascal triangle,
% Provide the kcd's in order 1-1 1-2 1-3 1-4...2-2 2-3 2-4...
% So we can be sloppy and just iterate through the kcd's in the loop below
if (length(kcd_num) == num_sets + pascalp(num_sets))
    kcd_vec = kcd_num;
elseif (length(kcd_num) == 1)
    kcd_vec = ones(1,(num_sets + pascalp(num_sets))).*kcd_num;
else
    error('compete_cyc: kcd should have length 1 or %d',num_sets + pascalp(num_sets));
end
kc_multi = max(kcd_vec);
if (length(dfA_num) == num_sets)
	dfA_vec = dfA_num;
	dfB_vec = dfA_num;
elseif (length(dfA_num) == 1)
	dfA_vec = ones(1,num_sets).*dfA_num;
	dfB_vec = dfA_vec;
else
    error('compete_cyc: dfA should have length 1 or %d',num_sets);
end
if numvarargs == 3
    if (length(dfB_num) == num_sets)
        dfB_vec = dfB_num;
    elseif (length(dfB_num) == 1)
        dfB_vec = ones(1,num_sets).*dfB_num;
else
    error('compete_cyc: dfB should have length 1 or %d',num_sets);
    end
end
%dfA_vec
%dfB_vec
% Mvec is a vector of monomer concentrations of length 3*num_sets
% MDD is a vector of dead-dead monomers, of length num_sets
% M_LiveA_LiveB
MLL = (1-dfA_vec) .* (1-dfB_vec) .* M_input_vec;
% M_LiveA_DeadB
MLD = (1-dfA_vec) .* (dfB_vec) .* M_input_vec;
% M_DeadA_LiveB
MDL = (dfA_vec) .* (1-dfB_vec) .* M_input_vec;
% Don't need to track dead-dead since concentration is constant
MDD = (dfA_vec) .* (dfB_vec) .* M_input_vec;
% Then reassort the monomers by class again

for isn = 1:num_sets
	Mvec((isn-1)*3+1) = MLL(isn);
	Mvec((isn-1)*3+2) = MLD(isn);
	Mvec((isn-1)*3+3) = MDL(isn);
end

% Build the matrices that we will use to calculate the delta concs in the
% integration routine
% There will be three monomers per set that can do anything 
num_monos = length(Mvec);
% The number of combinations seems like it might be
% = 4 * (num_monos + (num_monos-1) + ... + num_monos-(num_monos-1))) 
% = 4 * 21 for 6 monomers = 84 but many are not formed because the center
% ends must be alive to ligate.
% Number of combos within a set that can actually form is 3 unique self-ligation
% products for a molecule L-i-L with 2 good ends, 2 unique products for L-i-L
% ligated to either L-i-D for a total 4, and three for combining the two L-i-D
% pairs, giving 10 for each set.
% Then inter-class there are 4 unique L-i-L+L-j-L dimers, 2 each = 4
% L-i-L+L-j-D, 2 each = 4 L-j-L+L-i-D, 4 L-i-D L-j-D combos, summing to 16
% for each inter-class set.
% So for two classes there shouls be 10+16+10 = 36 combinations

% For three classes there should be 10+10+10+16+16+16 = 78 formable dimer
% combinations, for four there will be 4*10 + 6*16 = 136 etc.
% pascalp already defined above
%pascalp = [0 1 3 6 10 15 21 28 36 45 55 66 78 91 ]; % etc. inter-pairings are triangular numbers from the Pascal triangle,
% i.e. binomial coefficients, generated as nchoosek(n+1,2) if needed where n+1 is the second element
% Plus there will be one monomer circle per set
% Plus two i-i dimer circles per set and two i-j dimer circles for each 
% i != j pairing
% Plus a large number of multimers that we lump together for now. Could
% split out trimers, tetramers etc. 
%
% Species column vector will be monomers then dimers then monomer circles
% then dimer circles then multimers, and finally multimer circles
%
% num_species = num_sets*( 3 monomers + 10 intra-class dimers + 1 circle + 2 dimer circle)
%   + pascal triangle number*(16 inter-class dimer + 2 interclass DC) +
%   three kinds of "Multimers" (L-L, L-D, D-D)
%num_species = num_sets*3 + num_sets*(10+1+2) + pascalp(num_sets)*(16+2) + 3 % = 53 for num_sets = 2
num_bimols = num_sets*(10) + pascalp(num_sets)*(16);
num_unicirc = num_sets;
% There are ABAB and ABBA double circles and each can be formed from
% interclass pairs
num_doublec = num_sets*2 + pascalp(num_sets)*2; % But we will at least assume that double circles within a size class have the same kcd
% The multimers are split into LL, LD, and DD clumps
num_multi = 4;
num_species = num_monos + num_unicirc + num_doublec + num_bimols + num_multi; % check whether it adds up.
% Conc_blank allows us to have an element of the concentration matrix that
% is always zero, for convenience so we don't have to check logic every
% time on whether an element exists (i.e. as opposed to having a "0" in an
% index column). Ensure that conc_vec(conc_blank) is set to and remains 0.
conc_blank = num_species + 1; 
%
%conc_vec = zeros(num_species+1,1);
conc0_vec = zeros(num_species+1,1);
conc0_vec(1:num_monos) = Mvec;
% The parameter matrix is constructed to describe the formation of each individual species.
% Each row describes one species:
% monomers, then bimolecalar products, then uni circles, then double
% circles, then 3 kinds of multimers.
% At this point we can also set up the cyclization of all double-ended
% molecules, since there is only one partner in a cyclization. 
%
% The nine columns are:
% 1: Species number that converts to a unimolecular circle. Number of nonzero 
%    entries is 2 * the number of LL monomers, conc_blank otherwise
% 2: unimolecular cyclization rate constant for species indexed in column 1; + in the row
%    for the circle being formed, - in the row for the linear monomer
% 3: Species number #1 that converts to a double circle. There is just one
%    for an ABAB double circle.
% 4: Species number #2 that converts to a double circle, there are two for
%    ABBA and BAAB bimolecular products that give the same circle. Set this
%    to conc_blank for an ABAB circle
% 5: cyclization rate constant for the DC species; + in the row
%    for the circle being formed, - in the row for the linear monomer, assumed
%    the same for ABAB and ABBA/BAAB. If there were bending this would not be
%    true.
% 6: Species number #1 for the monomer that is incorporated into the bimolecular product for
%    the row, conc_blank otherwise
% 7: The other species that is in the bimolecular product
% 8: The number of ways the bimol product can be formed (= 2 for LA_isn_LB +
%    LA_isn_LB -> LA_isn_LB_LA_isn_LB, = 1 otherwise).
% 9: The number of live ends for each species -- this does not change, no
%    need to have a separate vector

param_matrix = zeros(num_species+1,9);
%live_ends_vec = zeros(num_species+1,1);

% Set up the parameter matrix by iterating through monomer sets 
% First do intra-set and cyclization, then pairwise interactions
% Each set of 3 monomers is LL LD DL (not including DD, which is stuck at a constant
%  concentration determined above); only LL can yield circles
for mono_num = 1:num_monos
    set_num = floor((mono_num+2)/3); % 111 222 333 444
    param_matrix(mono_num,9) = 1 + ~mod(mono_num+2,3); % 211 211 211 211
    conc0_vec(mono_num) = Mvec(mono_num);
% This enters kc in the appropriate cyclization formation matrix
% Should not need to iterate through circles separately
    if ~mod(mono_num+2,3) % 100 100 100 100
        % There is no way to form a monomer, and the decay is set up
        % separately 
        param_matrix(mono_num,:) = [ mono_num -1*kc_vec(set_num) conc_blank conc_blank 0 conc_blank conc_blank 0 2];
        % This does set up the complete parameter matrix row for the unimolecualr circle
        param_matrix(num_monos+num_bimols+set_num,:) = [ mono_num kc_vec(set_num) conc_blank conc_blank 0 conc_blank conc_blank 0 0];
% Next four lines are redundant here; written for parallel structure and
% completeness
%        conc0_vec(num_monos+num_bimols+set_num) = 0;
%        live_ends_vec(num_monos+num_bimols+set_num) = 0;
%    else
%        param_matrix(mono_num,:) = [ 0 0 0 0 0 0 0];
    else
        % Enter conc_blank in appropriate places so deltaconc program can
        % be written more simply
        param_matrix(mono_num,:) = [ conc_blank 0 conc_blank conc_blank 0 conc_blank conc_blank 0 1];
    end
end
% Any monomer can react with any other species, we take care of dead ends
% being dead later
%    breakdown_matrix(left_mono,:) = kb;
% Add the kc term for loss of the monomer for the monomer with two live ends
%    breakdown_matrix(left_mono,left_mono) = breakdown_matrix(left_mono,left_mono) + ~mod(left_mono+2,3)*kc_vec(set_num);
% Then iterate through the other monomers to set up the formation terms for individual bimolecular products
% list intra-class dimers LA-isn-LB_LA-isn-LB etc.
%  LA_1_LB_LA_1_LB, LA_1_LB_LA_1_DB, DA_1_LB_LA_1_LB, DA_1_LB_LA_1_DB,
%  LA_1_LB_LB_1_LA, LA_1_LB_LB_1_DA, DA_1_LB_LB_1_DA,
%  LB_1_LA_LA_1_LB, LB_1_LA_LA_1_DB, DB_1_LA_LA_1_DB

% List inter-class dimers

%  LA_1_LB_LA_2_LB, DA_1_LB_LA_2_LB, LA_1_LB_LA_2_DB, DA_1_LB_LA_2_DB, ends up as rows 1-4
%  LA_2_LB_LA_1_LB, DA_2_LB_LA_1_LB, LA_2_LB_LA_1_DB, DA_2_LB_LA_1_DB,
%  LA_1_LB_LB_2_LA, DA_1_LB_LB_2_LA, LA_1_LB_LB_2_DA, DA_1_LB_LB_2_DA,
%  LB_1_LA_LA_2_LB, DB_1_LA_LA_2_LB, LB_1_LA_LA_2_DB, DB_1_LA_LA_2_DB,

% We step through the bimolecular products set by set "triangle-wise" ie.
% 1+1 1+2 1+3 1+4 ... then 2+2 2+3 2+4 ... then 3+3 3+4 ...
bimol_base = num_monos;
kcd_index = 0;
dc_base = num_monos + num_bimols + num_sets;
for isn0 = 1:num_sets
    bi_inter_abab{isn0} = [];
end
for isn1 = 1:num_sets
% intra-class
% just increment bimol_num to keep track instead of doing it smarter
    if isn1 > 1
        bimol_base = bimol_base + 10 + (num_sets-isn1+1)*16;
        dc_base = dc_base + (num_sets-isn1+2)*2;
    end
    kcd_index = kcd_index + 1;
    bimol_num = bimol_base + 1; %1
% First the ABAB dimers
% The first intraclass bimol is LA-isn1-LB_LA_isn1_LB so it can cyclize to an ABAB dimer
% circle, which is the only way to get that circle. The bimol  is formed from a single monomer in two separate ways
    param_matrix(bimol_num,:) = [conc_blank 0 bimol_num conc_blank -kcd_vec(kcd_index) 3*(isn1-1)+1 3*(isn1-1)+1 2 2];
    bi_intra_abab{isn1} = [bimol_num:bimol_num+3];
    param_matrix(dc_base+1,:) = [ conc_blank 0 bimol_num conc_blank kcd_vec(kcd_index) conc_blank conc_blank 0 0];
    dc_abab(isn1) = dc_base+1;
%    live_ends_vec(num_monos+num_bimols+num_sets+1) = 0;
%    conc0_vec(num_monos+num_bimols+num_sets+1) = 0;
%	live_ends_vec(num_monos + bimol_num) = 2;  % Could calculate this by summing the values for the monomers - 2
% The second and third intraclass bimols have one dead end each and cannot
% cyclize
    bimol_num = bimol_num+1; %2  
    param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+1 3*(isn1-1)+2 1 1];
%	live_ends_vec(num_monos + bimol_num) = 1;
    bimol_num = bimol_num+1; %3
    param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+1 3*(isn1-1)+3 1 1];
%	live_ends_vec(num_monos + bimol_num) = 1;
% The fourth intraclass bimol has two dead ends, so it connot react further
    bimol_num = bimol_num+1; %4
    param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+2 3*(isn1-1)+3 1 0];
%	live_ends_vec(num_monos + bimol_num) = 0;
% And the ABBA
% First the LA-isn1-LB_LB_isn1_LA dimer
    bimol_num = bimol_num+1; %5
    dimer_circ_input1 = bimol_num;
    param_matrix(bimol_num,:) = [conc_blank 0 dimer_circ_input1 conc_blank -kcd_vec(kcd_index) 3*(isn1-1)+1 3*(isn1-1)+1 1 2];
%	live_ends_vec(num_monos + bimol_num) = 2;
    %param_matrix(num_monos+num_bimols+num_sets+(isn1-1)*2+2,:) = [ conc_blank 0 (num_monos + bimol_num) kcd_vec(kcd_index) 0 0 0 ];
    bimol_num = bimol_num+1; %6
    param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+1 3*(isn1-1)+3 1 1];
%	live_ends_vec(num_monos + bimol_num) = 1;
    bimol_num = bimol_num+1; %7
    param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+3 3*(isn1-1)+3 1 0];
%	live_ends_vec(num_monos + bimol_num) = 0;
% And the BAAB
    bimol_num = bimol_num+1; %8
    param_matrix(bimol_num,:) = [conc_blank 0 conc_blank bimol_num -kcd_vec(kcd_index) 3*(isn1-1)+1 3*(isn1-1)+1 1 2];
%	live_ends_vec(num_monos + bimol_num) = 2;
%   Now we can set up the second dimer circle
    dc_num = dc_base+2;
    param_matrix(dc_num,:) = [ conc_blank 0 dimer_circ_input1 bimol_num kcd_vec(kcd_index) conc_blank conc_blank 0 0 ];
    bimol_num = bimol_num+1; %9
    param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+1 3*(isn1-1)+2 1 1];
%	live_ends_vec(num_monos + bimol_num) = 1;
    bimol_num = bimol_num+1; %10
    param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+2 3*(isn1-1)+2 1 0];
%	live_ends_vec(num_monos + bimol_num) = 0;        
% inter-class. Without loss of generality we list the ABAB products in that orientation, and
% specify the lower-numbered class on the left for ABBA and BAAB
% List inter-class dimers
% The four entries in the first column will all form at the same rate but
% in theory the DC's could form differently so we keep track separately.

%  LA_1_LB_LA_2_LB, DA_1_LB_LA_2_LB, LA_1_LB_LA_2_DB, DA_1_LB_LA_2_DB, ends up as rows 1-4
%  LA_2_LB_LA_1_LB, DA_2_LB_LA_1_LB, LA_2_LB_LA_1_DB, DA_2_LB_LA_1_DB,
%  LA_1_LB_LB_2_LA, DA_1_LB_LB_2_LA, LA_1_LB_LB_2_DA, DA_1_LB_LB_2_DA,
%  LB_1_LA_LA_2_LB, DB_1_LA_LA_2_LB, LB_1_LA_LA_2_DB, DB_1_LA_LA_2_DB,

    for isn2 = isn1+1:num_sets
        kcd_index = kcd_index + 1; % So these must be provided in the proper order
% First the ABAB dimers
% First four have lower isn to the left
        bi_inter_abab{isn1} = [ bi_inter_abab{isn1} [bimol_num+5:bimol_num+8]];
        bi_inter_abab{isn2} = [ bi_inter_abab{isn2} [bimol_num+1:bimol_num+4]];
        bimol_num = bimol_num+1; %1
        param_matrix(bimol_num,:) = [conc_blank 0 bimol_num conc_blank -kcd_vec(kcd_index) 3*(isn1-1)+1 3*(isn2-1)+1 1 2];
        dimer_circ_input1 = bimol_num;
        bimol_num = bimol_num+1; %2
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+3 3*(isn2-1)+1 1 1];
        bimol_num = bimol_num+1; %3
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+1 3*(isn2-1)+2 1 1];
        bimol_num = bimol_num+1; %4
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+3 3*(isn2-1)+2 1 1];
%       second set of four have lower number on the right. They are different.
        bimol_num = bimol_num+1; %5
        param_matrix(bimol_num,:) = [conc_blank 0 bimol_num conc_blank -kcd_vec(kcd_index) 3*(isn2-1)+1 3*(isn1-1)+1 1 2];
        dc_num = dc_num + 1;
        param_matrix(dc_num,:) = [ conc_blank 0 dimer_circ_input1 bimol_num kcd_vec(kcd_index) conc_blank conc_blank 0 0];
        bimol_num = bimol_num+1; %6
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn2-1)+3 3*(isn1-1)+1 1 1];
        bimol_num = bimol_num+1; %7
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn2-1)+1 3*(isn1-1)+2 1 1];
        bimol_num = bimol_num+1; %8
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn2-1)+3 3*(isn1-1)+2 1 1];
% Then the ABBA dimers. Constructed so the lower isn is always to the left
        bimol_num = bimol_num+1; %9
        param_matrix(bimol_num,:) = [conc_blank 0 bimol_num conc_blank -kcd_vec(kcd_index) 3*(isn1-1)+1 3*(isn2-1)+1 1 2];
        dimer_circ_input1 = bimol_num;
        bimol_num = bimol_num+1; %10
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+3 3*(isn2-1)+1 1 1];
        bimol_num = bimol_num+1; %11
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+1 3*(isn2-1)+3 1 1];
        bimol_num = bimol_num+1; %12
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+3 3*(isn2-1)+3 1 1];
% Then the BAAB dimers. Constructed so the lower isn is always to the left
        bimol_num = bimol_num+1; %13
        param_matrix(bimol_num,:) = [conc_blank 0 bimol_num conc_blank -kcd_vec(kcd_index) 3*(isn1-1)+1 3*(isn2-1)+1 1 2];
        dc_num = dc_num + 1;
        param_matrix(dc_num,:) = [ conc_blank 0 dimer_circ_input1 bimol_num kcd_vec(kcd_index) conc_blank conc_blank 0 0];
        bimol_num = bimol_num+1; %14
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+2 3*(isn2-1)+1 1 1];
        bimol_num = bimol_num+1; %15
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+1 3*(isn2-1)+2 1 1];
        bimol_num = bimol_num+1; %16
        param_matrix(bimol_num,:) = [conc_blank 0 conc_blank conc_blank 0 3*(isn1-1)+2 3*(isn2-1)+2 1 1];

% bimol_form_instraset_mat = ...
%     [ 2 0 0 ; ... % First four lines are the ABAB dimers
%       1 1 0 ; ...
%       1 0 1 ; ...
%       0 1 1 ; ...
%       1 0 0 ; ... % Then the ABBA
%       1 0 1 ; ...
%       0 0 1 ; ...
%       1 0 0 ; ... % Then the BAAB
%       1 1 0 ; ...
%       0 1 0 ]
% bimol_form_interset_mat_half =  ...
%     [ 1 0 0 ; ... % First four lines are the ABAB dimers
%       0 0 1 ; ...
%       1 0 0 ; ...
%       0 0 1 ; ...
%       1 0 0 ; ... % rows 5-8 switches 1 and 2
%       1 0 0 ; ...
%       0 1 0 ; ...
%       0 1 0 ; ...
%       1 0 0 ; ... % ABBA dimers
%       0 0 1 ;
%       1 0 0 ;
%       
    end
end
      

      
% There is no mechanism for making monomers so the top num_monos row remain zeros        
% Then any pair of monomers can make a dimer        

% For example, for two sets the concs vector should contain the six types of monomer above that we
% track, plus
%  LA_1_LB,DA_1_LB,LA_1_DB,LA_2_LB,DA_1_LB,LA_1_DB
%
%  LA_1_LB_LA_1_LB, DA_1_LB_LA_1_LB, LA_1_LB_LA_1_DB, DA_1_LB_LA_1_DB,
%  LA_2_LB_LA_2_LB, DA_2_LB_LA_2_LB, LA_2_LB_LA_2_DB, DA_2_LB_LA_2_DB,
%  LA_1_LB_LA_2_LB, DA_1_LB_LA_2_LB, LA_1_LB_LA_2_DB, DA_1_LB_LA_2_DB,
%  LA_2_LB_LA_1_LB, DA_2_LB_LA_1_LB, LA_2_LB_LA_1_DB, DA_2_LB_LA_1_DB,
%
%  LA_1_LB_LB_1_LA, DA_1_LB_LB_1_LA, DA_1_LB_LB_1_DA,
%  LA_2_LB_LB_2_LA, DA_2_LB_LB_2_LA, DA_2_LB_LB_2_DA,
%  LA_1_LB_LB_2_LA, DA_1_LB_LB_2_LA, LA_1_LB_LB_2_DA, DA_1_LB_LB_2_DA,
%
%  LB_1_LA_LA_1_LB, DB_1_LA_LA_1_LB, DB_1_LA_LA_1_DB,
%  LB_2_LA_LA_2_LB, DB_2_LA_LA_2_LB, DB_2_LA_LA_2_DB,
%  LB_1_LA_LA_2_LB, DB_1_LA_LA_2_LB, LB_1_LA_LA_2_DB, DB_1_LA_LA_2_DB,
%  C1, C2,
%  DC_11_ABAB, DC_11_ABBA, DC_12_ABAB, DC_12_AABB, DC_22_ABAB, DC_22_AABB, 
%  and lump all multimers


% For debugging send the param matrix to the screen
%param_matrix;

% 12dimer, 11DC, 22DC, 12DC,multimers
%concs0 = [M1_ActiveA_ActiveB,M1_ActiveA_DeadB,M1_DeadA_ActiveB,M2_ActiveA_ActiveB,M2_ActiveA_DeadB,M2_DeadA_ActiveB,...
%    zeros(1,45)];
% Can use ode23 or ode15s. ode15s seems to be about 40% faster.
% Inputs to deltaconcs that are not in the param_matrix are...kb,num_monos,
% num_bimols, num_unicirc, num_doublec
num_vec = [num_monos num_bimols num_unicirc num_doublec num_multi num_species]
% Try to pre-calculate everything that does not change
monos_indices = 1:num_vec(1);
%param_matrix(monos_indices,1)
%conc0_vec(param_matrix(monos_indices,1))
%conc0_vec(param_matrix(monos_indices,1)).*param_matrix(monos_indices,2)
bb = (num_vec(1)+1):(num_vec(1)+num_vec(2));
%param_matrix(bb,1)
unic_indices = (num_vec(1)+num_vec(2)+1):(num_vec(1)+num_vec(2)+num_vec(3));
%param_matrix(unic_indices,1)
dc_indices = (num_vec(1)+num_vec(2)+num_vec(3)+1):(num_vec(1)+num_vec(2)+num_vec(3)+num_vec(4));
%param_matrix(dc_indices,1)
multi_indices = (num_vec(1)+num_vec(2)+num_vec(3)+num_vec(4))+1:(num_vec(1)+num_vec(2)+num_vec(3)+num_vec(4)+num_vec(5)) ;
%param_matrix(multi_indices,1)
% Multimers can come from monomer + bimol or bimol + bimol or
% multimer+monomer or multimer + bimol or multimer + multimer
% How to keep track? This can be doen in the dconcs routine by just adding
% everything.
% We do need to separately have a line for cyclization of the double-live
% ended multimer to give the fourth multimer population

param_matrix(multi_indices(1),:) = [multi_indices(1) -kc_multi conc_blank conc_blank 0 conc_blank conc_blank 1 2];
param_matrix(multi_indices(4),:) = [multi_indices(1) kc_multi conc_blank conc_blank 0 conc_blank conc_blank 1 0];
        
% For the multi-calc we need to have lists of monomers with 2 active ends,
% bimols with 2 ends, monos with 1 end, bimols with 1 end.
% I think this clumsy logic is correct now.
[temp_ind_2,~] = find(param_matrix(:,9) == 2) ;
[temp_ind_1,~] = find(param_matrix(:,9) == 1) ;
monos_two_active = temp_ind_2(find(temp_ind_2 <= num_monos)) ;
bimols_two_active = temp_ind_2(find(temp_ind_2 > num_monos & temp_ind_2 < (num_monos + num_bimols))) ;
monos_one_active = temp_ind_1(find(temp_ind_1 <= num_monos)) ;
bimols_one_active = temp_ind_1(find(temp_ind_1 > num_monos & temp_ind_1 < (num_monos + num_bimols))) ;

%bimols_two_active = param_matrix(bb,9) == 2 
%bimols_one_active = param_matrix(bb,9) == 1 
%index_set = [monos_indices,bb,unic_indices,dc_indices,multi_indices,monos_two_active,monos_one_active,bimols_two_active,bimols_one_active]
% test delta concs
%concs = conc0_vec
%testconcs(3,concs,param_matrix,...
%    monos_indices,bb,unic_indices,dc_indices,multi_indices,monos_two_active,monos_one_active,bimols_two_active,bimols_one_active,num_species,kb)
options = odeset('RelTol',1e-7,'AbsTol',1e-9,'NonNegative',1:num_species);
%
[t,concs] = ode15s(@(t,concs) deltaconcs(t,concs,param_matrix,...
    monos_indices,bb,unic_indices,dc_indices,multi_indices,monos_two_active,monos_one_active,bimols_two_active,bimols_one_active,num_species,kb),...
    [0 end_time],conc0_vec',options);
%concs
%products = [M0-concs(:,3)-concs(:,4) concs(:,3) concs(:,4)];
% For plotting we want to see the monomers, the unimolecular circles, the intra-class
% ABAB dimers that end up indistinguishable from circles, and the dimer circles that end up indistinguishable
figure
hold on
% Still lumping the multi concentrations even though some of them are
% circular
conc_linear_multis = sum(concs(:,multi_indices(1:3)),2);
conc_circ_multis = sum(concs(:,multi_indices(4)),2);
for kplot = 1:num_sets
    conc_monomers(:,kplot) = sum(concs(:,[3*kplot-2:3*kplot]),2) + MDD(kplot);
    conc_bimol_intra_ABAB(:,kplot) = sum(concs(:,bi_intra_abab{kplot}),2);
    conc_bimol_inter_ABAB(:,kplot) = sum(concs(:,bi_inter_abab{kplot}),2);
    conc_unicircle(:,kplot) = concs(:,unic_indices(kplot));
    conc_dc_ABAB(:,kplot) = 2*concs(:,dc_abab(kplot));
%    all_else (kplot) = 
    plot(t,conc_monomers(:,kplot),'Color',[.25+.5*kplot/num_sets .25+.5*kplot/num_sets .25+.5*kplot/num_sets],'LineWidth',1); % gray
    plot(t,conc_bimol_intra_ABAB(:,kplot),'Color',[.25+.75*kplot/num_sets 0 0],'LineWidth',1); % red
    plot(t,conc_bimol_inter_ABAB(:,kplot),'Color',[.25+.75*kplot/num_sets 0 .25+.75*kplot/num_sets],'LineWidth',1); % red
    plot(t,conc_unicircle(:,kplot),'Color',[0 .25+.75*kplot/num_sets 0],'LineWidth',2); % green
    plot(t,conc_dc_ABAB(:,kplot),'Color',[0 0 .25+.75*kplot/num_sets],'LineWidth',5);% blue   
end
plot(t,conc_linear_multis,'Color',[.25+.75*kplot/num_sets 0 .25+.75*kplot/num_sets],'LineWidth',2);% pink
plot(t,conc_circ_multis,'Color',[ 0 .25+.75*kplot/num_sets .25+.75*kplot/num_sets],'LineWidth',2);% teal
ax = axis;
box on

% The second plot just shows everything
figure
hold on
plot(t,concs(:,monos_indices),'Color',[0.5 0.5 0.5]); % gray
plot(t,concs(:,bb),'Color',[1 0 0]); % red
plot(t,concs(:,unic_indices),'Color',[0 1 0]); % green
plot(t,concs(:,dc_indices),'Color',[0 0 1]); % blue
plot(t,concs(:,multi_indices),'Color',[0.75 0.25 0]); % brown
legend
axis(ax);
%products = concs;
end

function dconcs = deltaconcs(~,concs,param_matrix,...
    monos_indices,bb,unic_indices,dc_indices,multi_indices,monos_two_active,monos_one_active,bimols_two_active,bimols_one_active,num_species,kb)
% here is where all the cleverness with the param_matrix comes into play
%
% the dot product below is essentially the concentration of all reactive
% ends in the mixture multiplied by kb. It is a scalar.
react_ends_rate = kb*concs'*(param_matrix(:,9));
M_2_conc_sum = sum(concs(monos_two_active));
B_2_conc_sum = sum(concs(bimols_two_active));
M_1_conc_sum = sum(concs(monos_one_active));
B_1_conc_sum = sum(concs(bimols_one_active));
dconcs = zeros(num_species+1,1);
%for ii = monos_indices -> looks like it can be done with an array
%statement.
% We assume that the param matrix entries have already identified which
% monomers make uni circles so we don't have to check the live ends
% number in the first line
dconcs(monos_indices) = concs(param_matrix(monos_indices,1)).*param_matrix(monos_indices,2)...
        - (concs(monos_indices).*param_matrix(monos_indices,9))*react_ends_rate ...
        - kb*concs(monos_indices).*concs(monos_indices).*param_matrix(monos_indices,9).^2;
dconcs(bb) = (concs(param_matrix(bb,3))+concs(param_matrix(bb,4))).*param_matrix(bb,5) ...
	+ kb*param_matrix(bb,8).*concs(param_matrix(bb,6)).*concs(param_matrix(bb,7)) ...
	- (concs(bb).*param_matrix(bb,9))*react_ends_rate ...
    - kb*concs(bb).*concs(bb).*param_matrix(bb,9).^2 ;
dconcs(unic_indices) = concs(param_matrix(unic_indices,1)).*param_matrix(unic_indices,2);
dconcs(dc_indices) = (concs(param_matrix(dc_indices,3))+concs(param_matrix(dc_indices,4))).*param_matrix(dc_indices,5);
% Third line below -- when two double-live multis ligate together, the
% product is one double-ended multi so the concentration decreases
% The last line represents cyclization of multimer
dconcs(multi_indices(1)) = 4*kb*(M_2_conc_sum + 0.5*B_2_conc_sum)*B_2_conc_sum +  ...
                            + 4*kb*(M_2_conc_sum + B_2_conc_sum)*concs(multi_indices(1)) ...
                            - 4*kb*concs(multi_indices(1))*concs(multi_indices(1)) ...
                            - 2*kb*(M_1_conc_sum + B_1_conc_sum)*concs(multi_indices(1)) ...
                            - 2*kb*concs(multi_indices(1))*concs(multi_indices(2)) ...
                            + concs(param_matrix(multi_indices(1),1))*param_matrix(multi_indices(1),2);
dconcs(multi_indices(2)) = 2*kb*(M_2_conc_sum + B_2_conc_sum)*B_1_conc_sum +  ...
                            + 2*kb*(M_2_conc_sum + B_2_conc_sum)*concs(multi_indices(2)) ...
                            + 2*kb*(M_1_conc_sum + B_1_conc_sum)*concs(multi_indices(1)) ...
                            + 2*kb*concs(multi_indices(1))*concs(multi_indices(2)) ...
                            - kb*(M_1_conc_sum + B_1_conc_sum)*concs(multi_indices(2)) ...
                            - kb*concs(multi_indices(2))*concs(multi_indices(2));                        
dconcs(multi_indices(3)) = kb*(M_1_conc_sum + 0.5*B_1_conc_sum)*B_1_conc_sum +  ...
                            + kb*(M_1_conc_sum + B_1_conc_sum)*concs(multi_indices(2)) ...
                            + kb*concs(multi_indices(2))*concs(multi_indices(2));
dconcs(multi_indices(4)) = concs(param_matrix(multi_indices(4),1))*param_matrix(multi_indices(4),2);                      
%sum(dconcs)
end

% function dconcs = testconcs(~,concs,param_matrix,...
%     monos_indices,bb,unic_indices,dc_indices,multi_indices,monos_two_active,monos_one_active,bimols_two_active,bimols_one_active,num_species,kb)
% % here is where all the cleverness with the param_matrix comes into play
% %
% % the dot product below is essentially the concentration of all reactive
% % ends in the mixture multiplied by kb. It is a scalar.
% react_ends_rate = kb*concs'*(param_matrix(:,9));
% M_2_conc_sum = sum(concs(monos_two_active));
% B_2_conc_sum = sum(concs(bimols_two_active));
% M_1_conc_sum = sum(concs(monos_one_active));
% B_1_conc_sum = sum(concs(bimols_one_active));
% dconcs = zeros(num_species+1,1);
% %for ii = monos_indices -> looks like it can be done with an array
% %statement.
% % We assume that the param matrix entries have already identified which
% % monomers make uni circles so we don't have to check the live ends
% % number in the first line
% dconcs(monos_indices) = concs(param_matrix(monos_indices,1)).*param_matrix(monos_indices,2);
% dconcs(monos_indices) = dconcs(monos_indices) - concs(monos_indices)*react_ends_rate ;
% dconcs(monos_indices) = dconcs(monos_indices)- kb*concs(monos_indices).*concs(monos_indices).*param_matrix(monos_indices,9).^2;
% dconcs(bb) = (concs(param_matrix(bb,3))+concs(param_matrix(bb,4))).*param_matrix(bb,5) ...
% 	+ kb*param_matrix(bb,8).*concs(param_matrix(bb,6)).*concs(param_matrix(bb,7)) ...
% 	- (concs(bb).*param_matrix(bb,9))*react_ends_rate ...
%     - kb*concs(bb).*concs(bb).*param_matrix(bb,9).^2 ;
% dconcs(unic_indices) = concs(param_matrix(unic_indices,1)).*param_matrix(unic_indices,2);
% dconcs(dc_indices) = (concs(param_matrix(dc_indices,3))+concs(param_matrix(dc_indices,4))).*param_matrix(dc_indices,5);
% % Third line below -- when two double-live multis ligate together, the
% % product is one double-ended multi so the concentration decreases
% dconcs(multi_indices(1)) = 4*kb*(M_2_conc_sum + 0.5*B_2_conc_sum)*B_2_conc_sum +  ...
%                             + 4*kb*(M_2_conc_sum + B_2_conc_sum)*concs(multi_indices(1)) ...
%                             - 4*kb*concs(multi_indices(1))*concs(multi_indices(1)) ...
%                             - 2*kb*(M_1_conc_sum + B_1_conc_sum)*concs(multi_indices(1)) ...
%                             - 2*kb*concs(multi_indices(1))*concs(multi_indices(2));
% dconcs(multi_indices(2)) = 2*kb*(M_2_conc_sum + B_2_conc_sum)*B_1_conc_sum +  ...
%                             + 2*kb*(M_2_conc_sum + B_2_conc_sum)*concs(multi_indices(2)) ...
%                             + 2*kb*(M_1_conc_sum + B_1_conc_sum)*concs(multi_indices(1)) ...
%                             + 2*kb*concs(multi_indices(1))*concs(multi_indices(2)) ...
%                             - kb*(M_1_conc_sum + B_1_conc_sum)*concs(multi_indices(2)) ...
%                             - kb*concs(multi_indices(2))*concs(multi_indices(2));                        
% dconcs(multi_indices(3)) = kb*(M_1_conc_sum + 0.5*B_1_conc_sum)*B_1_conc_sum +  ...
%                             + kb*(M_1_conc_sum + B_1_conc_sum)*concs(multi_indices(2)) ...
%                             + kb*concs(multi_indices(2))*concs(multi_indices(2));
% end
% 

