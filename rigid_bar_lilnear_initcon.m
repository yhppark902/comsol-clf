clear all
clc
% Model Import
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model')
model.modelPath(['D:\usr']);

%% FEM parameters

model.param.set('M0', '0[N*m]');
model.param.set('dtFem', '0.02[s]');
model.param.set('phi', '0[rad]');
model.param.set('phi_t', '0[rad/s]');

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', [-0.025 0]);
model.component('comp1').geom('geom1').feature('r1').set('size', [0.05 1]);
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').run('pt1');
model.component('comp1').geom('geom1').create('rot1', 'Rotate');
model.component('comp1').geom('geom1').feature('rot1').selection('input').set({'r1'});
model.component('comp1').geom('geom1').feature('rot1').set('rot', 180);
model.component('comp1').geom('geom1').run;

model.component('comp1').physics.create('mbd', 'MultibodyDynamics', 'geom1');
model.component('comp1').physics('mbd').create('rd1', 'RigidDomain', 2);
model.component('comp1').physics('mbd').feature('rd1').selection.set([1]);
model.component('comp1').physics('mbd').feature('rd1').feature('crp1').selection.set([3]);
model.component('comp1').physics('mbd').create('hgj1', 'HingeJoint', -1);
model.component('comp1').physics('mbd').feature('hgj1').set('Source', 'fixed');
model.component('comp1').physics('mbd').feature('hgj1').set('Destination', 'rd1');
model.component('comp1').physics('mbd').feature('hgj1').create('afm1', 'AppliedForceAndMoment', -1);
model.component('comp1').physics('mbd').create('gr1', 'Gravity', 2);
model.component('comp1').physics('mbd').feature('gr1').selection.set([1]);

model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').probe.create('dom1', 'Domain');
model.component('comp1').probe.create('dom2', 'Domain');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').label('Structural steel');
model.component('comp1').material('mat1').propertyGroup('def').set('density', '7850[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', '200e9[Pa]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', '0.30');

model.component('comp1').physics('mbd').prop('PhysicsSymbols').set('PhysicsSymbols', true);
model.component('comp1').physics('mbd').prop('d').set('d', 0.05);
model.component('comp1').physics('mbd').feature('rd1').set('CenterOfRotationType', 'CentroidOfSelectedEntities');
model.component('comp1').physics('mbd').feature('rd1').set('EntityLevel', 'Point');
model.component('comp1').physics('mbd').feature('rd1').set('ConsistentInitialization', 'ForceInitialValues');
model.component('comp1').physics('mbd').feature('hgj1').set('CenterOfJointType', 'UserDefined');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('AppliedOnSelectedAttachment', 'Joint');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('PointOfApplicationType', 'CenterOfJoint');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('Mz', 'input1');
model.component('comp1').physics('mbd').feature('hgj1').feature('afm1').set('Ms', 'M0');
model.component('comp1').physics('mbd').feature('rd1').set('InitialValueType', 'locallyDefined');
model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phi', 'phi');
model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phit', 'phi_t');

model.component('comp1').mesh('mesh1').run;
model.component('comp1').probe('dom1').set('expr', 'mbd.rd1.phi');
model.component('comp1').probe('dom2').set('intsurface', true);
model.component('comp1').probe('dom2').set('intvolume', true);
model.component('comp1').probe('dom2').set('expr', 'mbd.rd1.th_tz');

model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.result.dataset.create('dset2', 'Solution');
model.result.dataset.create('avh1', 'Average');
model.result.dataset('dset2').set('probetag', 'dom1');
model.result.dataset('avh1').set('probetag', 'dom1');
model.result.dataset('avh1').set('data', 'dset2');
model.result.dataset('avh1').selection.geom('geom1', 2);
model.result.dataset('avh1').selection.set([1]);
model.result.numerical.create('pev1', 'EvalPoint');
model.result.numerical('pev1').set('probetag', 'dom1');
model.result.create('pg1', 'PlotGroup1D');
model.result('pg1').create('tblp1', 'Table');
model.result('pg1').feature('tblp1').set('probetag', 'dom1');
model.component('comp1').probe('dom1').genResult([]);

model.study('std1').feature('time').set('tlist', 'range(0,dtFem/3,dtFem)');


model.sol('sol1').attach('std1');

%% run
model.sol('sol1').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('control', 'user');
model.sol('sol1').feature('v1').set('initsol', 'zero');
model.sol('sol1').feature('v1').set('initmethod', 'init');
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('v1').feature('comp1_u').set('scalemethod', 'auto');
model.sol('sol1').feature('v1').feature('comp1_u').set('resscalemethod', 'auto');
model.sol('sol1').feature('v1').feature('comp1_mbd_rd1_u').set('scalemethod', 'auto');
model.sol('sol1').feature('v1').feature('comp1_mbd_rd1_phi').set('scalemethod', 'auto');
model.sol('sol1').feature('t1').label('Time-Dependent Solver 1.1');
model.sol('sol1').feature('t1').set('control', 'user');
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'onevery');
model.sol('sol1').runAll;

model.result.numerical('pev1').set('const', {'mbd.refpntx' '0' 'Reference point for moment computation, x component'; 'mbd.refpnty' '0' 'Reference point for moment computation, y component'; 'mbd.refpntz' '0' 'Reference point for moment computation, z component'});
model.result.numerical('pev1').setResult;
model.result('pg1').set('xlabel', 'Time (s)');
model.result('pg1').set('ylabel', 'Rigid body rotation (rad), Domain Probe 1');
model.result('pg1').feature('tblp1').label('Probe Table Graph 1');
model.result('pg1').feature('tblp1').set('plotcolumninput', 'manual');
model.result('pg1').feature('tblp1').set('legend', true);
model.result('pg1').feature('tblp1').set('table', 'tbl1');
model.result('pg1').feature('tblp1').set('plotcolumns', [2]);
model.result('pg1').run;
%% System Identification
t=mphglobal(model,'root.t');
array_t=t;
th=mphglobal(model,'mbd.hgj1.th');
array_th=th;
th_t=mphglobal(model,'mbd.hgj1.th_t');
array_th_t=th_t;
array_t_ode23=[];
array_th_ode23=[];


% model.sol('sol1').feature('v1').set('initmethod', 'sol');
% model.sol('sol1').feature('v1').set('initsol', 'sol1');
% model.sol('sol1').feature('v1').set('solnum', 'last');
M = mphstate(model,'sol1','out',{'Mc' 'MA' 'MB' 'A' 'B' 'C' 'D' 'x0', 'Null', 'ud'},'input', {'M0'}, 'output', {'comp1.dom1','comp1.dom2'}, 'sparse', 'off', 'initmethod','init');

%% Assumed feedback gian, CLF
kp =6;
kd =5;
clf_rate=3;
slack = 1e3;
u_max=7;
u_ref=pi;
nt=2;
phi=0;
phi_t=0;
x_init=[0;0];

for k=1:10
func = @(tt,xx)M.MA*xx+M.MB*100;
opt=odeset('mass',M.Mc,'jacobian',M.MA);
[t_ode,x_ode]=ode23s(func,[0 0.02],x_init,opt);
y=M.C*x_ode';
x_init=x_ode(end,:);

array_t_ode23=[array_t_ode23;0.02*(1+k)+t_ode];
array_th_ode23=[array_th_ode23;y'];

model.param.set('M0', '100[N*m]');
t=mphglobal(model,'root.t')+t(end);
array_t=[array_t;t];
th=mphglobal(model,'mbd.hgj1.th');
array_th=[array_th;th];
phi=array_th(end);
th_t=mphglobal(model,'mbd.hgj1.th_t');
array_th_t=[array_th_t;th_t];
phi_t=array_th_t(end);
model.param.set('phi', strcat( num2str(array_th(end)),'[rad]'));
model.param.set('phi_t', strcat( num2str(array_th_t(end)),'[rad/s]'));
model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phi', 'phi');
model.component('comp1').physics('mbd').feature('rd1').feature('init1').set('phit', 'phi_t');
model.sol('sol1').runAll;
end

figure(1)
plot(array_t,array_th)


hold on

plot(array_t_ode23,array_th_ode23,'--')
legend('comsol','ode23s');
xlabel('time(s)');
ylabel('theat(rad)');