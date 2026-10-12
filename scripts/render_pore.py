"""Run with Blender's --python, then -- --data ... --mesh ... --style ... --out ..."""
from pathlib import Path
import argparse,json,sys
import numpy as np
import bpy
from mathutils import Vector,Matrix

def linear(h):
    x=[int(h[i:i+2],16)/255 for i in [1,3,5]]
    return tuple(v/12.92 if v<.04045 else ((v+.055)/1.055)**2.4 for v in x)

def cut_surface(V,F,origin,normal):
    d=(V-origin)@normal;vertices=V.tolist();faces=F[(d[F]<=0).all(1)].tolist();cache={}
    def intersection(i,j):
        if abs(d[i])<1e-10:return int(i)
        if abs(d[j])<1e-10:return int(j)
        key=tuple(sorted((int(i),int(j))))
        if key not in cache:cache[key]=len(vertices);vertices.append((V[i]+d[i]/(d[i]-d[j])*(V[j]-V[i])).tolist())
        return cache[key]
    for triangle in F[(d[F]<=0).any(1)&(d[F]>0).any(1)]:
        polygon=[]
        for k in range(3):
            i,j=triangle[k],triangle[(k+1)%3]
            if d[i]<=0:polygon.append(int(i))
            if (d[i]<=0)!=(d[j]<=0):polygon.append(intersection(i,j))
        for k in range(1,len(polygon)-1):
            if len({polygon[0],polygon[k],polygon[k+1]})==3:faces.append([polygon[0],polygon[k],polygon[k+1]])
    F=np.array(faces,np.int32);used,inverse=np.unique(F,return_inverse=True)
    return np.array(vertices)[used],inverse.reshape(-1,3)

def collection(name):
    c=bpy.data.collections.new(name);bpy.context.scene.collection.children.link(c);return c

def material(name,h,roughness):
    m=bpy.data.materials.new(name);m.use_nodes=True;b=m.node_tree.nodes.get('Principled BSDF')
    b.inputs['Base Color'].default_value=(*linear(h),1);b.inputs['Roughness'].default_value=roughness;b.inputs['Coat Weight'].default_value=.1
    return m

def sphere_mesh(name,m,resolution):
    bpy.ops.mesh.primitive_uv_sphere_add(segments=resolution,ring_count=resolution//2,radius=1)
    o=bpy.context.object;mesh=o.data;mesh.name=name;mesh.materials.append(m)
    for p in mesh.polygons:p.use_smooth=True
    bpy.data.objects.remove(o,do_unlink=True);return mesh

def cylinder_mesh(name,m):
    bpy.ops.mesh.primitive_cylinder_add(vertices=12,radius=1,depth=1)
    o=bpy.context.object;mesh=o.data;mesh.name=name;mesh.materials.append(m)
    for p in mesh.polygons:p.use_smooth=len(p.vertices)==4
    bpy.data.objects.remove(o,do_unlink=True);return mesh

def stick(name,p,q,r,mesh,col):
    p,q=Vector(p),Vector(q);d=q-p;o=bpy.data.objects.new(name,mesh);col.objects.link(o)
    o.location=(p+q)/2;o.scale=(r,r,d.length);o.rotation_mode='QUATERNION';o.rotation_quaternion=d.to_track_quat('Z','Y')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--data',required=True);parser.add_argument('--mesh',required=True);parser.add_argument('--style',required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--mode',choices=['full','cutaway'],default='full')
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    data=json.loads(Path(args.data).read_text());style=json.loads(Path(args.style).read_text());z=np.load(args.mesh)
    if data['pore_components']!=1:raise ValueError('Display has not been validated as one pore')
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    bpy.ops.wm.read_factory_settings(use_empty=True);scene=bpy.context.scene
    frame=collection('01 Framework');guest=collection('02 Guests');pores=collection('03 Pore');studio=collection('04 Studio')
    atoms=data['framework_atoms'];X=np.array([a['position_A'] for a in atoms]);V=z['vertices_A'];F=z['triangles']
    d=np.array(style['camera_direction'],float);d/=np.linalg.norm(d)
    axis=np.array(data['pore_axis_uvw'])@np.array(data['lattice_A']);axis/=np.linalg.norm(axis)
    # 把侧视方向投影到孔轴的法平面；不同 CIF 的笛卡尔基底不会改变侧视含义。
    d-=axis*np.dot(d,axis)
    if np.linalg.norm(d)<1e-8:raise ValueError('Requested side camera is parallel to pore axis; choose another direction')
    d/=np.linalg.norm(d)
    if args.mode=='cutaway':
        if abs(np.dot(d,axis))>1e-5:raise ValueError('Longitudinal cut requires camera direction perpendicular to verified pore axis')
        origin=np.array(data['guest_anchor_position_A']);V,F=cut_surface(V,F,origin,d)
        np.savez_compressed(out/'cutaway_mesh.npz',vertices_A=V,triangles=F,plane_origin_A=origin,plane_normal=d)
    elements=set(a['element'] for a in atoms);guest_elements=set(a['element'] for g in data['guests'] for a in g['atoms'])
    if not elements|guest_elements<=style['palette'].keys():raise ValueError('Add explicit element colors to the style; no random palette substitutions')
    if not elements<=style['framework_radii_A'].keys():raise ValueError('Add explicit framework display radii')
    mats={e:material('Framework '+e,style['palette'][e],.46) for e in elements}
    meshes={e:sphere_mesh('Framework sphere '+e,mats[e],16) for e in elements}
    cylinders={e:cylinder_mesh('Framework bond '+e,mats[e]) for e in elements}
    # 框架统一细球棍；所有位置从显示数据读取，不强化特定配体。
    for i,a in enumerate(atoms):
        e=a['element'];o=bpy.data.objects.new(f'Framework {i} {e}',meshes[e]);frame.objects.link(o)
        o.location=a['position_A'];o.scale=(style['framework_radii_A'][e],)*3;o['source_position_A']=a['position_A'];o['source_index']=a['source_index'];o['framework_atom']=True
    for k,(i,j,_) in enumerate(data['framework_bonds']):
        p,q=X[i],X[j];mid=(p+q)/2;r=style.get('framework_bond_radius_A',.033)*(.67 if 'H' in [atoms[i]['element'],atoms[j]['element']] else 1)
        for n,start,end,e in [(0,p,mid,atoms[i]['element']),(1,mid,q,atoms[j]['element'])]:stick(f'Framework bond {k}.{n}',start,end,r,cylinders[e],frame)
    gm={e:material('Guest '+e,style.get('guest_palette',{}).get(e,style['palette'][e]),.30) for e in guest_elements}
    gs={e:sphere_mesh('Guest sphere '+e,gm[e],48) for e in guest_elements};gc={e:cylinder_mesh('Guest bond '+e,gm[e]) for e in guest_elements}
    guest_bounds=[]
    # 客体保留真实坐标；spacefill scale 只改变显示半径，不改变计算半径。
    for i,g in enumerate(data['guests']):
        for j,a in enumerate(g['atoms']):
            e=a['element'];r=(data['vdw_radii_A'][e]*float(style.get('guest_vdw_scale',.55))) if style.get('guest_representation')=='spacefill' else style['guest_ball_radii_A'][e]
            o=bpy.data.objects.new(f'Guest {i}.{j} {e}',gs[e]);guest.objects.link(o);o.location=a['position_A'];o.scale=(r,)*3
            o['source_position_A']=a['position_A'];o['source_index']=a['source_index'];o['guest_atom']=True;o['guest_molecule']=i;o['display_radius_A']=r
            p=np.array(a['position_A']);guest_bounds.extend([p-r,p+r])
        if style.get('guest_representation')!='spacefill':
            for k,(i1,i2,_) in enumerate(g['bonds']):
                a,b=g['atoms'][i1],g['atoms'][i2];p=np.array(a['position_A']);q=np.array(b['position_A']);mid=(p+q)/2;r=style.get('guest_bond_radius_A',.078)
                stick(f'Guest {i} bond {k}.0',p,mid,r,gc[a['element']],guest);stick(f'Guest {i} bond {k}.1',mid,q,r,gc[b['element']],guest)
    mesh=bpy.data.meshes.new('Contact Surface');mesh.from_pydata(V.tolist(),[],F.tolist());mesh.update()
    for p in mesh.polygons:p.use_smooth=True
    pore=bpy.data.objects.new('Single computed pore',mesh);pores.objects.link(pore);pore.visible_shadow=False;pore['probe_radius_A']=data['probe_radius_A']
    mat=bpy.data.materials.new('Pore inside outside transparency');mat.use_nodes=True;nt=mat.node_tree;nt.nodes.clear()
    geom=nt.nodes.new('ShaderNodeNewGeometry');color=nt.nodes.new('ShaderNodeMixRGB')
    # 内外层分别用同色系浅青蓝与蓝灰；透明着色避免折射导致的暗原子。
    color.inputs[1].default_value=(*linear(style['pore_outer']),1);color.inputs[2].default_value=(*linear(style['pore_inner']),1);nt.links.new(geom.outputs['Backfacing'],color.inputs[0])
    b=nt.nodes.new('ShaderNodeBsdfPrincipled');nt.links.new(color.outputs[0],b.inputs['Base Color']);b.inputs['Roughness'].default_value=.43;b.inputs['Coat Weight'].default_value=.1
    transparent=nt.nodes.new('ShaderNodeBsdfTransparent');mix=nt.nodes.new('ShaderNodeMixShader');mix.name='Pore opacity';opacity=style['cutaway_opacity'] if args.mode=='cutaway' else style['pore_opacity'];mix.inputs[0].default_value=opacity
    nt.links.new(transparent.outputs[0],mix.inputs[1]);nt.links.new(b.outputs[0],mix.inputs[2]);output=nt.nodes.new('ShaderNodeOutputMaterial');nt.links.new(mix.outputs[0],output.inputs['Surface']);mesh.materials.append(mat)
    if style.get('section_outline',False):raise ValueError('This renderer preserves the approved no-outline style; an outline requires a deliberate separate implementation.')
    requested_up=np.array(style['camera_up'],float);u=np.cross(d,axis)
    if np.dot(u,requested_up)<0:u=-u
    u/=np.linalg.norm(u);right=np.cross(u,d)
    bounds=np.vstack([X,V,*([np.array(guest_bounds)] if guest_bounds else [])]);center=(bounds.min(0)+bounds.max(0))/2;P=bounds-center
    width,height=style.get('resolution',[2000,1300]);aspect=width/height
    # Blender 正交尺度按横向宽度拟合，兼顾高度，给完整原子和客体留白。
    scale=max(np.ptp(P@right)+2,(np.ptp(P@u)+2)*aspect)*1.07
    cd=bpy.data.cameras.new('Side orthographic');cam=bpy.data.objects.new('Side orthographic',cd);studio.objects.link(cam)
    cd.type='ORTHO';cd.sensor_fit='HORIZONTAL';cd.ortho_scale=scale;cd.clip_end=1000;cam.location=Vector(center+d*90);cam.rotation_mode='QUATERNION';cam.rotation_quaternion=Matrix((right.tolist(),u.tolist(),d.tolist())).transposed().to_quaternion();scene.camera=cam
    world=bpy.data.worlds.new('Bright studio');scene.world=world;world.use_nodes=True;world.node_tree.nodes['Background'].inputs['Strength'].default_value=.55
    world.node_tree.nodes['Background'].inputs['Color'].default_value=(1,1,1,1);span=float(np.linalg.norm(bounds.max(0)-bounds.min(0)))
    # 宽主光、填充光与柔和轮廓光，使孔道有立体感，原子明亮。
    for name,offset,energy,size in [('Key',[-.7,.9,1.2],5.,.6),('Front fill',[.85,.1,1],2.3,1.6),('Left fill',[-1,-.1,.45],1.5,1.75),('Lower fill',[0,-1,.55],1.2,2),('Rim',[.45,.85,-1],3.5,.9)]:
        ld=bpy.data.lights.new(name,'AREA');ld.energy=energy*span**2;ld.shape='DISK';ld.size=size*span;light=bpy.data.objects.new(name,ld);studio.objects.link(light)
        light.location=Vector(center+(right*offset[0]+u*offset[1]+d*offset[2])*span*.72);light.rotation_euler=(Vector(center)-light.location).to_track_quat('-Z','Y').to_euler()
    scene.render.engine='CYCLES';scene.cycles.samples=style.get('samples',128);scene.cycles.use_denoising=True;scene.cycles.transparent_max_bounces=32;scene.cycles.adaptive_threshold=.015
    backend='CPU';prefs=bpy.context.preferences.addons['cycles'].preferences
    for kind in (['METAL'] if sys.platform=='darwin' else ['OPTIX','CUDA','HIP']):
        try:
            prefs.compute_device_type=kind;prefs.get_devices();devices=[device for device in prefs.devices if device.type==kind]
            if devices:
                for device in prefs.devices:device.use=device.type==kind
                scene.cycles.device='GPU';backend=kind;break
        except (TypeError,RuntimeError):pass
    if backend=='CPU':scene.cycles.device='CPU'
    scene.render.film_transparent=True;scene.render.resolution_x=width;scene.render.resolution_y=height;scene.render.resolution_percentage=100;scene.render.image_settings.file_format='PNG';scene.render.image_settings.color_mode='RGBA'
    scene.view_settings.view_transform='Standard';scene.view_settings.look='None';scene['probe_radius_A']=data['probe_radius_A'];scene['pore_method']='Periodic Contact Surface; see calculation.json for native network validation';scene['view']=args.mode
    scene.render.filepath=str(out/'figure_transparent.png');path=out/'figure.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(path));bpy.ops.wm.open_mainfile(filepath=str(path));bpy.ops.render.render(write_still=True)
    scene=bpy.context.scene;objects=[o for o in scene.objects if o.get('framework_atom') or o.get('guest_atom')]
    error=max(np.linalg.norm(np.array(o.location)-np.array(o['source_position_A'])) for o in objects)
    pore=bpy.data.objects['Single computed pore'];mesh_error=np.max(np.linalg.norm(np.array([v.co[:] for v in pore.data.vertices])-V,axis=1))
    if error>1e-4 or mesh_error>1e-4:raise ValueError('Saved scene positions differ from validated input')
    if any(o.type=='FONT' or o.type=='CURVE' or 'BFDH' in o.name for o in scene.objects):raise ValueError('Unexpected text/outline/morphology object')
    report={'passed':True,'saved_scene_reopened':True,'rendered_from_saved_scene':True,'framework_atoms':len(atoms),'guest_molecules':len(data['guests']),'guest_atoms':sum(bool(o.get('guest_atom')) for o in scene.objects),'position_error_A':float(error),'mesh_error_A':float(mesh_error),'pore_triangles':len(F),'mode':args.mode,'backend':backend,'opacity':opacity,'section_outline':False,'guest_representation':style.get('guest_representation'),'guest_vdw_scale':style.get('guest_vdw_scale')}
    (out/'render_validation.json').write_text(json.dumps(report,indent=2)+'\n');(out/'style_used.json').write_text(json.dumps(style,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
