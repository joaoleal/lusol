function [methodinfo,structs,enuminfo,ThunkLibName]=libclusol_proto
%LIBCLUSOL_PROTO Create structures to define interfaces found in 'clusol'.

%This function was generated by loadlibrary.m parser version 1.1.6.35 on Mon Oct 21 10:56:54 2013
%perl options:'clusol.i -outfile=libclusol_proto.m -thunkfile=libclusol_thunk_glnxa64.c -header=clusol.h'
ival={cell(1,0)}; % change 0 to the actual number of functions to preallocate the data.
structs=[];enuminfo=[];fcnNum=1;
fcns=struct('name',ival,'calltype',ival,'LHS',ival,'RHS',ival,'alias',ival,'thunkname', ival);
MfilePath=fileparts(mfilename('fullpath'));
ThunkLibName=fullfile(MfilePath,'libclusol_thunk_glnxa64');
% void clu1fac ( int64_t * m , int64_t * n , int64_t * nelem , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * iploc , int64_t * iqloc , int64_t * ipinv , int64_t * iqinv , double * w , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu1fac'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu6sol ( int64_t * mode , int64_t * m , int64_t * n , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu6sol'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu8rpc ( int64_t * mode1 , int64_t * mode2 , int64_t * m , int64_t * n , int64_t * jrep , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform , double * diag , double * vnorm ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8rpc'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr'};fcnNum=fcnNum+1;
% void clu6mul ( int64_t * mode , int64_t * m , int64_t * n , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu6mul'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu8adc ( int64_t * mode , int64_t * m , int64_t * n , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform , double * diag , double * vnorm ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8adc'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr'};fcnNum=fcnNum+1;
% void clu8adr ( int64_t * m , int64_t * n , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform , double * diag ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8adr'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr'};fcnNum=fcnNum+1;
% void clu8dlc ( int64_t * m , int64_t * n , int64_t * jdel , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8dlc'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu8dlr ( int64_t * mode , int64_t * m , int64_t * n , int64_t * idel , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8dlr'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu8mod ( int64_t * mode , int64_t * m , int64_t * n , double * beta , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8mod'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu8rpr ( int64_t * mode1 , int64_t * mode2 , int64_t * m , int64_t * n , int64_t * irep , double * v , double * w , double * wnew , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8rpr'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
methodinfo=fcns;