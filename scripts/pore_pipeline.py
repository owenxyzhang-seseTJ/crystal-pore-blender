"""Ordered CIF -> periodic Contact Surface -> single-pore display geometry.

Coordinates are always fractional @ the same row-stored display lattice.
No electronic-structure optimization, occupancy selection or guest relocation.
"""
from pathlib import Path
import argparse,itertools,json,subprocess
import numpy as np
from scipy.ndimage import label
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from pymatgen.core import Structure,Lattice
from periodic_connectivity import periodic_components

VDW={'H':1.09,'C':1.70,'N':1.55,'O':1.52,'S':1.80,'Fe':2.,'Co':2.,'Zn':1.39}
BONDS={('C','H'):1.20,('C','C'):1.78,('C','N'):1.68,('C','O'):1.56,('C','S'):2.02,('C','Fe'):2.25,('H','N'):1.20,('H','O'):1.20,('H','S'):1.45,('C','F'):1.50,('C','Cl'):1.92,('C','Br'):2.05,('C','P'):1.95,('O','P'):1.75}

def dump(p,d):
    Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')

def load_job(path):
    path=Path(path).resolve();j=json.loads(path.read_text())
    if 'cif' in j:j['cif']=str((path.parent/j['cif']).resolve())
    return j

def intake(cif,A=None):
    source=Structure.from_file(cif)
    if any(not s.is_ordered or abs(sum(s.species.values())-1)>1e-6 for s in source):
        raise ValueError('Disorder/partial occupancy: supply an explicitly ordered model; source unchanged.')
    lattice=np.array(source.lattice.matrix if A is None else A,float)
    if lattice.shape!=(3,3) or abs(np.linalg.det(lattice))<1e-8:raise ValueError('Invalid display lattice')
    if not np.allclose(lattice@lattice.T,source.lattice.matrix@source.lattice.matrix.T,atol=1e-4,rtol=1e-6):
        raise ValueError('Display lattice changes cell metrics; only consistent Cartesian reorientation is allowed here.')
    atoms=[{'source_index':i,'element':s.specie.symbol,'label':s.label,'fractional':(s.frac_coords%1).tolist()} for i,s in enumerate(source)]
    return source,{'source_cif':str(Path(cif).resolve()),'lattice_A':lattice.tolist(),'atoms':atoms,'ordered':True}

def radii(job,elements):
    r={**VDW,**job.get('vdw_radii_A',{})}
    missing=set(elements)-r.keys()
    if missing:raise ValueError(f'Explicit vdW radii required for {sorted(missing)}')
    if any(float(r[e])<=0 for e in elements):raise ValueError('vdW radii must be positive')
    return r

def replicas(A,reach):
    # 用倒易基底得到晶胞面间距，确保强倾斜晶胞也覆盖所需邻近球。
    heights=1/np.linalg.norm(np.linalg.inv(A),axis=0)
    pads=np.maximum(1,np.ceil(reach/heights).astype(int))
    return np.array(list(itertools.product(*[range(-n,n+1) for n in pads])))

def fields(F,elements,A,r,probe,step,max_voxels):
    N=np.ceil(np.linalg.norm(A,axis=1)/step).astype(int)
    if np.prod(N)>max_voxels:raise ValueError(f'Grid {N.tolist()} exceeds max_voxels; choose an appropriate cell/resolution.')
    frac=np.indices(N).reshape(3,-1).T/N;X=frac@A
    margin=max(np.linalg.norm(A,axis=1)/N)*1.8
    clearance=np.full(len(X),probe+margin)
    for e in sorted(set(elements)):
        atomF=F[np.array(elements)==e];reach=r[e]+probe+margin;T=replicas(A,reach)
        P=(atomF[None,:,:]+T[:,None,:]).reshape(-1,3)@A
        distance=cKDTree(P).query(X,distance_upper_bound=reach)[0]
        clearance=np.minimum(clearance,distance-r[e])
    clearance=clearance.reshape(N).astype(np.float32)
    _,connected,components=periodic_components(clearance>=probe)
    if not connected.any():raise ValueError('No percolating accessible-center set for this probe/radii model.')
    centers=frac[connected.ravel()];T=replicas(A,probe+margin)
    P=(centers[None,:,:]+T[:,None,:]).reshape(-1,3)@A
    distance=cKDTree(P).query(X,distance_upper_bound=probe+margin)[0]
    distance=np.minimum(distance,probe+margin)
    contact=(probe-distance).reshape(N).astype(np.float32)
    return clearance,connected,contact,components

def compute(job,out):
    source,data=intake(job['cif'],job.get('lattice_A'));A=np.array(data['lattice_A'])
    groups=job.get('guest_groups',[]);removed=[i for group in groups for i in group]
    if len(removed)!=len(set(removed)) or any(i<0 or i>=len(source) for i in removed):raise ValueError('Invalid/overlapping guest indices')
    retained=[i for i in range(len(source)) if i not in removed]
    if not retained:raise ValueError('No framework remains')
    F=np.array([data['atoms'][i]['fractional'] for i in retained]);els=[data['atoms'][i]['element'] for i in retained]
    r=radii(job,els);probe=float(job.get('probe_radius_A',.9));step=float(job.get('grid_step_A',.14))
    if probe<=0 or step<=0:raise ValueError('Positive probe and grid step required')
    _,mask,contact,components=fields(F,els,A,r,probe,step,int(job.get('max_voxels',12000000)))
    out.mkdir(parents=True,exist_ok=True);dump(out/'intake.json',data)
    framework=Structure(source.lattice,[source[i].specie for i in retained],[source[i].frac_coords for i in retained])
    framework.to(filename=str(out/'framework_for_pores.cif'))
    np.savez_compressed(out/'periodic_probe_field.npz',accessible_centers=mask,contact_level_A=contact,lattice_A=A,probe_radius_A=probe,framework_fractional=F,framework_elements=np.array(els))
    report={'probe_radius_A':probe,'grid_counts':list(map(int,contact.shape)),'max_grid_step_A':float(max(np.linalg.norm(A,axis=1)/np.array(contact.shape))),'vdw_radii_A':{e:r[e] for e in sorted(set(els))},'source_atoms':len(source),'framework_atoms':len(retained),'removed_guest_indices':removed,'hydrated_or_source_lattice_preserved':True,'surface_definition':'Contact Surface','contact_fraction':float((contact>0).mean()),'periodic_components':components,'clearance_note':'Accessible-center decisions exact within required cutoff; clearance outside the cutoff is capped, not usable for maximum cavity diameter.','zeopp':{'ran':False}}
    if 'convergence_step_A' in job:
        _,_,coarse,_=fields(F,els,A,r,probe,float(job['convergence_step_A']),int(job.get('max_voxels',12000000)))
        delta=abs(float((coarse>0).mean())-report['contact_fraction'])
        report['coarse_contact_fraction']=float((coarse>0).mean());report['fraction_resolution_difference']=delta
        if delta>float(job.get('convergence_tolerance',.025)):raise ValueError('Contact volume has not converged; refine grid before drawing.')
    executable=job.get('zeopp_executable')
    if executable:
        if any(abs(r[e]-VDW.get(e,float('nan')))>1e-8 or e not in VDW for e in set(els)):
            raise ValueError('Native Zeo++ radii must match overrides; configure matching native radii separately before claiming agreement.')
        command=[str(Path(executable).resolve()),'-ha','-res',str((out/'zeopp.res').resolve()),'-chan',str(probe),str((out/'zeopp.chan').resolve()),str((out/'framework_for_pores.cif').resolve())]
        result=subprocess.run(command,text=True,capture_output=True,timeout=300)
        (out/'zeopp.log').write_text(result.stdout+'\n'+result.stderr)
        good=result.returncode==0 and (out/'zeopp.res').exists() and (out/'zeopp.chan').exists()
        report['zeopp']={'ran':True,'returncode':result.returncode,'passed':good,'command':command}
        dump(out/'calculation.json',report)
        if not good:raise ValueError('Native Zeo++ calculation failed; inspect zeopp.log. No successful native result claimed.')
    report['passed']=True;dump(out/'calculation.json',report)
    print(json.dumps({'calculation':str(out),'contact_fraction':report['contact_fraction'],'components':components}))

def gentle_smooth(V,F):
    edges=np.sort(np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1)
    unique,count=np.unique(edges,axis=0,return_counts=True);boundary=np.unique(unique[count==1])
    matrix=coo_matrix((np.ones(len(unique)*2),(np.r_[unique[:,0],unique[:,1]],np.r_[unique[:,1],unique[:,0]])),shape=(len(V),len(V))).tocsr()
    degree=np.asarray(matrix.sum(1)).ravel();X=V.copy()
    # 三轮温和平滑，固定原有开放边界；结果仅供显示。
    for _ in range(3):
        for strength in [.25,-.26]:
            delta=matrix@X/degree[:,None]-X;delta[boundary]=0;X+=strength*delta
    return X,unique

def write_obj(path,V,F):
    with path.open('w') as f:
        for p in V:f.write('v '+' '.join(f'{v:.7f}' for v in p)+'\n')
        for triangle in F:f.write('f '+' '.join(str(int(i)+1) for i in triangle)+'\n')

def bond_limit(e1,e2,job):
    node=set(job['node_elements']);donors=set(job.get('donor_elements',['N','O']))
    if (e1 in node and e2 in donors) or (e2 in node and e1 in donors):return float(job.get('coordination_cutoff_A',2.60))
    key=tuple(sorted((e1,e2)))
    overrides={tuple(sorted(k.split('-'))):v for k,v in job.get('bond_cutoffs_A',{}).items()}
    return overrides.get(key,BONDS.get(key,0))

def templates(data,job):
    A=np.array(data['lattice_A']);lattice=Lattice(A);atoms=data['atoms'];n=len(atoms)
    F=np.array([a['fractional'] for a in atoms]);guest=set(i for group in job.get('guest_groups',[]) for i in group)
    adjacency={i:[] for i in range(n) if i not in guest and atoms[i]['element'] not in job['node_elements']}
    for i,j in itertools.combinations(adjacency,2):
        limit=bond_limit(atoms[i]['element'],atoms[j]['element'],job)
        if not limit:continue
        distance,image=lattice.get_distance_and_image(F[i],F[j])
        if .35<distance<=limit:
            adjacency[i].append((j,np.array(image,int)));adjacency[j].append((i,-np.array(image,int)))
    components=[];seen=set()
    for start in adjacency:
        if start in seen:continue
        shifts={start:np.zeros(3,int)};stack=[start]
        while stack:
            i=stack.pop();seen.add(i)
            for j,T in adjacency[i]:
                proposed=shifts[i]+T
                if j in shifts:
                    if not np.array_equal(shifts[j],proposed):raise ValueError('Infinite covalent component: whole-ligand crop is inapplicable. Use a periodic network display crop.')
                else:shifts[j]=proposed;stack.append(j)
        ids=sorted(shifts);components.append([(i,F[i]+shifts[i]) for i in ids])
    return components

def display(job,calculation,out):
    from skimage.measure import marching_cubes
    data=json.loads((calculation/'intake.json').read_text());report=json.loads((calculation/'calculation.json').read_text())
    if not report.get('passed'):raise ValueError('Calculation not validated')
    z=np.load(calculation/'periodic_probe_field.npz');A=z['lattice_A'];C=z['contact_level_A'];N=np.array(C.shape)
    _,current=intake(job['cif'],job.get('lattice_A'))
    if current['atoms']!=data['atoms'] or not np.allclose(current['lattice_A'],A,atol=1e-5):raise ValueError('CIF/coordinates do not match cached field intake')
    removed=sorted(i for g in job.get('guest_groups',[]) for i in g)
    if removed!=sorted(report['removed_guest_indices']):raise ValueError('Guest selection differs from cached calculation')
    expected_radii=radii(job,[a['element'] for i,a in enumerate(data['atoms']) if i not in removed])
    if any(abs(expected_radii[e]-v)>1e-8 for e,v in report['vdw_radii_A'].items()) or abs(float(job.get('probe_radius_A',.9))-float(z['probe_radius_A']))>1e-8:raise ValueError('Probe/radii do not match cached calculation')
    atoms=data['atoms'];lattice=Lattice(A);guests=[]
    for record in job.get('guest_images',[]):
        ids=job['guest_groups'][record['group']];anchor=ids[0];fa=np.array(atoms[anchor]['fractional']);shift=np.array(record['shift'],int)
        molecular=[]
        for i in ids:
            f=np.array(atoms[i]['fractional']);_,T=lattice.get_distance_and_image(fa,f);f=f+T+shift
            molecular.append({'source_index':i,'element':atoms[i]['element'],'position_A':(f@A).tolist(),'fractional':f.tolist()})
        guest_bonds=[]
        for i,j in itertools.combinations(range(len(molecular)),2):
            distance=np.linalg.norm(np.array(molecular[i]['position_A'])-np.array(molecular[j]['position_A']))
            if .35<distance<=bond_limit(molecular[i]['element'],molecular[j]['element'],job):guest_bonds.append([i,j,float(distance)])
        guests.append({'source_guest':record['group'],'periodic_translation':shift.tolist(),'atoms':molecular,'bonds':guest_bonds})
    if 'crop_fractional' not in job:raise ValueError('Explicit crop_fractional min/max is required to select one pore')
    lo=np.array(job['crop_fractional']['min']);hi=np.array(job['crop_fractional']['max'])
    if not np.all(hi>lo):raise ValueError('Invalid crop bounds')
    start=np.floor(lo*N).astype(int);stop=np.ceil(hi*N).astype(int);axes=[np.arange(a,b+1) for a,b in zip(start,stop)]
    field=C[np.ix_(*[x%N[i] for i,x in enumerate(axes)])].copy();labels,n=label(field>0)
    anchors=np.array([g['atoms'][0]['fractional'] for g in guests])
    if not len(anchors):anchors=np.array([job['pore_anchor_fractional']])
    indices=np.rint(anchors*N-start).astype(int)
    if np.any(indices<0) or np.any(indices>=np.array(field.shape)):raise ValueError('Guest/pore anchor outside crop')
    keep=[int(labels[tuple(i)]) for i in indices]
    if 0 in keep or len(set(keep))!=1:raise ValueError('Anchors are not inside ONE pore component; select joint periodic images/crop.')
    selected=labels==keep[0]
    end_axes=job.get('open_crop_axes',[0])
    for d in set(range(3))-set(end_axes):
        if any(np.take(selected,k,axis=d).any() for k in [0,-1]):raise ValueError('Crop truncates lateral pore wall; expand lateral bounds without selecting another pore.')
    field[(field>0)&~selected]=-.14
    V,F,_,_=marching_cubes(field,level=0,allow_degenerate=False);V=((V+start)/N)@A
    smooth,edges=gentle_smooth(V,F) if job.get('smooth_display',True) else (V.copy(),np.unique(np.sort(np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0))
    graph=coo_matrix((np.ones(len(edges)*2),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(V),len(V)))
    if connected_components(graph,directed=False,return_labels=False)!=1:raise ValueError('Cropped surface contains more than one connected mesh')
    motion=np.linalg.norm(smooth-V,axis=1)
    if motion.max()>float(job.get('max_display_smoothing_A',.10)):raise ValueError('Display smoothing exceeds allowed motion')
    r=radii(job,[a['element'] for a in atoms]);tree=cKDTree(V);records=[];components=templates(data,job)
    inv=np.linalg.inv(A);cutoff=float(job.get('wall_extra_A',.45))
    for group in components:
        ids=[i for i,f in group];base=np.array([f for i,f in group]);padding=(max(r.values())+cutoff)*np.linalg.norm(inv,axis=0)
        lower=np.floor(lo-base.max(0)-padding).astype(int);upper=np.ceil(hi-base.min(0)+padding).astype(int)
        for T in itertools.product(*[range(a,b+1) for a,b in zip(lower,upper)]):
            X=(base+T)@A;heavy=np.array([atoms[i]['element']!='H' for i in ids])
            if not heavy.any() or not np.any(tree.query(X[heavy])[0]<np.array([r[atoms[i]['element']]+cutoff for i in ids if atoms[i]['element']!='H'])):continue
            for i,f,x in zip(ids,base+T,X):records.append({'source_index':i,'element':atoms[i]['element'],'position_A':x.tolist(),'fractional':f.tolist(),'role':'ligand'})
    if not records:raise ValueError('No nearby framework ligands selected')
    donors=np.array([a['position_A'] for a in records if a['element'] in job.get('donor_elements',['N','O'])]);donor_tree=cKDTree(donors) if len(donors) else None
    for atom in atoms:
        if atom['element'] not in job['node_elements']:continue
        f=np.array(atom['fractional']);padding=3*np.linalg.norm(inv,axis=0)
        lower=np.floor(lo-f-padding).astype(int);upper=np.ceil(hi-f+padding).astype(int)
        for T in itertools.product(*[range(a,b+1) for a,b in zip(lower,upper)]):
            ff=f+T;x=ff@A
            if donor_tree is not None and donor_tree.query(x)[0]<=float(job.get('coordination_cutoff_A',2.60)):
                records.append({'source_index':atom['source_index'],'element':atom['element'],'position_A':x.tolist(),'fractional':ff.tolist(),'role':'node'})
    X=np.array([a['position_A'] for a in records]);els=[a['element'] for a in records]
    if cKDTree(X).query_pairs(.2):raise ValueError('Duplicate/overlapping displayed framework sites')
    bonds=[]
    maxcut=max([2.65,float(job.get('coordination_cutoff_A',2.60)),*job.get('bond_cutoffs_A',{}).values()])
    for i,j in sorted(cKDTree(X).query_pairs(maxcut)):
        d=np.linalg.norm(X[i]-X[j])
        if .35<d<=bond_limit(els[i],els[j],job):bonds.append([i,j,float(d)])
    coordinate_error=0.
    for atom in records+[a for g in guests for a in g['atoms']]:
        difference=np.array(atom['fractional'])-np.array(data['atoms'][atom['source_index']]['fractional'])
        coordinate_error=max(coordinate_error,float(np.linalg.norm((difference-np.rint(difference))@A)))
    if coordinate_error>1e-5:raise ValueError('Display coordinates have been changed beyond integer lattice translations')
    out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'pore_mesh.npz',raw_vertices_A=V,vertices_A=smooth,triangles=F)
    write_obj(out/'pore_raw.obj',V,F);write_obj(out/'pore_display.obj',smooth,F)
    axis=job.get('pore_axis_uvw')
    if axis is None:raise ValueError('Set pore_axis_uvw from verified periodic winding/native channel analysis')
    dump(out/'display.json',{'framework_atoms':records,'framework_bonds':bonds,'guests':guests,'lattice_A':A.tolist(),'pore_axis_uvw':axis,'probe_radius_A':float(z['probe_radius_A']),'vdw_radii_A':r,'guest_anchor_position_A':(anchors@A).mean(0).tolist(),'pore_components':1})
    dump(out/'validation.json',{'passed':True,'framework_atoms':len(records),'guest_molecules':len(guests),'one_pore_component':True,'anchor_labels':keep,'lateral_wall_not_clipped':True,'coordinate_modulo_lattice_error_A':coordinate_error,'smoothing_max_motion_A':float(motion.max()),'smoothing_mean_motion_A':float(motion.mean()),'triangles':len(F)})
    print(json.dumps({'display':str(out),'atoms':len(records),'guest_molecules':len(guests),'smoothing_max_A':float(motion.max())}))

def main():
    parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('inspect');p.add_argument('--cif',required=True);p.add_argument('--out',type=Path,required=True)
    for name in ['compute','display']:
        p=sub.add_parser(name);p.add_argument('--job',required=True);p.add_argument('--out',type=Path,required=True)
        if name=='display':p.add_argument('--calculation',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='inspect':
        _,data=intake(args.cif);args.out.mkdir(parents=True,exist_ok=True);dump(args.out/'intake.json',data);print('Indexed atoms:',len(data['atoms']))
    elif args.command=='compute':compute(load_job(args.job),args.out.resolve())
    else:display(load_job(args.job),args.calculation.resolve(),args.out.resolve())

if __name__=='__main__':main()
