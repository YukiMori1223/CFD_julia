import Plots
using Printf
using Base.Iterators
using SparseArrays
using Profile
using Plots
@enum BOUNDARY begin
    BC_NONE    =0
    BC_WALL    =1
    BC_INLET   =2
    BC_OUTLET  =3
end

function plot_velocity_staggered(u,v,w,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min)
    z_slice=2
    xval =round.(collect(x_min : dx :x_min+dx*(Nx)),digits=3)
    yval =round.(collect(y_min : dx :y_min+dy*(Nx)),digits=3)

    #グラフ表示
    scale=0.1
    x_rangeX=collect(x_min : dx :x_min+dx*(Nx))
    y_rangeX=collect(y_min-dy*0.5 : dy :y_min+dy*(Ny+0.5))
    x_fig, y_fig = (repeat(x_rangeX, outer=length(y_rangeX)), repeat(y_rangeX, inner=length(x_rangeX)))
    Zero_x=zeros(Nx+1, Ny+2, Nz+2)
    #表示する座標
    vel_fig=Plots.quiver(
                    x_fig,y_fig,quiver=(scale*((u))[:,:,z_slice][:],Zero_x[:]),
                    line_z=repeat(u[:,:,z_slice][:], inner=4),
                    color=Plots.cgrad(:redsblues),
                    xticks=(xval,xval),
                    yticks=(yval,yval),
                    title="velocity")

    x_rangeY=collect(x_min-dx*0.5: dx :x_min+dx*(Nx+0.5))
    y_rangeY=collect(y_min: dy :y_min+dy*(Ny))
    x_fig, y_fig = (repeat(x_rangeY, outer=length(y_rangeY)), repeat(y_rangeY, inner=length(x_rangeY)))
    Zero_y=zeros(Nx+2, Ny+1, Nz+2)
    #表示する座標
    vel_fig=Plots.quiver!(
                    x_fig,y_fig,quiver=(Zero_y[:],scale*((v))[:,:,z_slice][:]),
                    line_z=repeat(v[:,:,z_slice][:], inner=4),
                    color=Plots.cgrad(:redsblues),
                    xticks=(xval,xval),
                    yticks=(yval,yval),
                    title="velocity",aspect_ratio =1.0,size=(800,800))
    return vel_fig
end

function plot_velocity_mixture(u,v,w,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min)
    z_slice=1
    xval =round.(collect(x_min +dx*0.5: dx :x_min+dx*(Nx-0.5)),digits=3)
    yval =round.(collect(y_min +dx*0.5: dx :y_min+dy*(Ny-0.5)),digits=3)

    u_cent=zeros(Nx, Ny, Nz)
    v_cent=zeros(Nx, Ny, Nz)
    for i in 1:Nx,j in 1:Ny,k in 1:Nz
        u_cent[i,j,k]=(u[i,j+1,k+1]+u[i+1,j+1,k+1])*0.5
        v_cent[i,j,k]=(v[i+1,j,k+1]+v[i+1,j+1,k+1])*0.5
    end        
    #グラフ表示
    scale=0.1
    x_fig, y_fig = (repeat(xval, outer=length(yval)), repeat(yval, inner=length(xval)))
    #表示する座標
    xaxis =round.(collect(x_min : dx :x_min+dx*(Nx+1)),digits=3)
    yaxis =round.(collect(y_min : dx :y_min+dy*(Ny+1)),digits=3)
    vel_mag=sqrt.(u_cent.^2+v_cent.^2)
    vel_fig=Plots.quiver(
                    x_fig,y_fig,quiver=(scale*u_cent[:,:,z_slice][:],scale*v_cent[:,:,z_slice][:]),
                    line_z=repeat(vel_mag[:,:,z_slice][:], inner=4),
                    color=Plots.cgrad(:darkrainbow),
                    #xticks=(xaxis,xaxis),
                    #yticks=(yaxis,yaxis),
                    xlim=(x_min-0.01,x_min+dx*(Nx+1)),
                    ylim=(y_min-0.01,y_min+dy*(Ny+1)),
                    title="velocity",aspect_ratio =1.0,
                    #arrow=Plots.arrow(:open)
                    #marker =:circle
                    )
    return vel_fig
end

function plot_velocity_contour(u,v,w,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min,x_range,y_range)
    z_slice=1
    xval =round.(collect(x_min +dx*0.5: dx :x_min+dx*(Nx-0.5)),digits=3)
    yval =round.(collect(y_min +dx*0.5: dx :y_min+dy*(Ny-0.5)),digits=3)

    u_cent=zeros(Nx, Ny, Nz)
    v_cent=zeros(Nx, Ny, Nz)
    for i in 1:Nx,j in 1:Ny,k in 1:Nz
        u_cent[i,j,k]=(u[i,j+1,k+1]+u[i+1,j+1,k+1])*0.5
        v_cent[i,j,k]=(v[i+1,j,k+1]+v[i+1,j+1,k+1])*0.5
    end        
    #グラフ表示
    scale=0.4
    #表示する座標
    xaxis =round.(collect(x_min : dx :x_min+dx*(Nx+1)),digits=3)
    yaxis =round.(collect(y_min : dx :y_min+dy*(Ny+1)),digits=3)
    vel_mag=sqrt.(u_cent.^2+v_cent.^2)
    vel_fig = Plots.heatmap( x_range,y_range,vel_mag[:,:,1]',xlabel = "x", ylabel = "y",
                    #xticks=(xaxis,xaxis),
                    #yticks=(yaxis,yaxis),
                    colorbar=:bottom,
                    color=Plots.cgrad(:lightrainbow),
                    colorbar_title="velocity (m/s)",
                    clims=(0.0,1.0)
                    #xlim=(x_min-0.01,x_min+dx*(Nx+1)),
                    #ylim=(y_min-0.01,y_min+dy*(Ny+1)),
                    )
    return vel_fig
end

#output
function output(u,v,w,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min,pres) 
    
    #x,y方向それぞれのセル中心座標
    x_range=collect(x_min+dx*0.5 : dx :x_min+dx*(Nx-0.5))
    y_range=collect(y_min+dy*0.5 : dy :y_min+dy*(Ny-0.5))
    
    vel_fig=plot_velocity_contour(u,v,w,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min,x_range,y_range)
    xaxis =round.(collect(x_min : dx :x_min+dx*(Nx+1)),digits=3)
    yaxis =round.(collect(y_min : dx :y_min+dy*(Ny+1)),digits=3)
    pres_fig=Plots.contourf( x_range,y_range,pres[:,:,1]',xlabel = "x", ylabel = "y",colorbar_title="Presuure (Pa)")

    range=collect(x_min+dx*0.5 : dx :x_min+dx*(Nx-0.5))            

    exp_yy=[0.65928,0.57492,0.51117,0.46604,0.33304,0.18719,0.05702,-0.0608,-0.10648,-0.27805,-0.38289,-0.2973,-0.2222,-0.20196,-0.18109]
    exp_yx=[0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547]  
    exp_xx=[0.9688, 0.9609,0.9531,0.9453,0.9063,0.8594,0.8047,0.5,0.2344,0.2266,0.1563,0.0938,0.0781,0.0703,0.0625]
    exp_xy=[-0.21388,-0.27669,-0.33714,-0.39188,-0.5155,-0.42665,-0.31966,0.02526,0.32235,0.33075,0.37095,0.32627,0.30353,0.29012,0.274850]
    x_vel_fig=Plots.plot(range,u[100,2:Ny+1,2],label="simulation")
    x_vel_fig=Plots.scatter!(exp_yx,exp_yy,label="Ghia, et al.",msw=0.0, ms=6,legend=:topleft,xlabel = "y position (m)",ylabel = "velocity u (m/s)")
    y_vel_fig=Plots.plot(range,v[2:Ny+1,100,2],label="simulation")
    y_vel_fig=Plots.scatter!(exp_xx,exp_xy,label="Ghia, et al.",msw=0.0, ms=6,xlabel = "x position (m)",ylabel="velocity v (m/s)")
    
    blank=Plots.plot(title="",grid=false, showaxis=false,xticks=false,yticks=false)
    l=@layout [
        e{0.01h}
        e{0.03w} a e{0.02w} b  
        e{0.01h}
        e{0.03w} c e{0.02w} d 
        e{0.01h}
        ]
    
    return Plots.plot(  blank, blank,vel_fig    , blank, pres_fig, blank,
                        blank,x_vel_fig  , blank, y_vel_fig ,blank ,title=["" "" "velocity" "" "presure" "" "" "U in center" "" "V in center" ""],
                        size=(1800,1400), layout=l)
end

function bound_vel!(u_f,v_f,w_f,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min,BC_type,u_BC,v_BC,w_BC)
    #X
    for j in 1:Ny+2,k in 1:Nz+2
        #X-
        BC_xminus=BC_type[1,j,k]
        #垂直方向
        if BC_xminus==BC_INLET
            u_f[1,j,k]=u_BC[1,j,k]
        elseif  BC_xminus==BC_OUTLET || BC_xminus==BC_NONE
            u_f[1,j,k]=u_f[2,j,k]
        elseif  BC_xminus==BC_WALL
            u_f[1,j,k]=0.0
        end
        #水平方向
        if j!=Ny+2
            if  BC_xminus==BC_OUTLET || BC_xminus==BC_NONE
                v_f[1,j,k]=v_f[2,j,k]
            elseif  BC_xminus==BC_WALL
                v_f[1,j,k]=2*v_BC[1,j,k]-v_f[2,j,k]
            end
        end     
        if k!=Nz+2
            if  BC_xminus==BC_OUTLET || BC_xminus==BC_NONE
                w_f[1,j,k]=w_f[2,j,k]
            elseif  BC_xminus==BC_WALL
                w_f[1,j,k]=2*w_BC[1,j,k]-w_f[2,j,k]
            end
        end

        #X+
        BC_xplus=BC_type[Nx+2,j,k]
        #垂直方向
        if BC_xplus==BC_INLET
            u_f[Nx+1,j,k]=u_BC[Nx+1,j,k]
        elseif  BC_xplus==BC_OUTLET || BC_xplus==BC_NONE
            u_f[Nx+1,j,k]=u_f[Nx,j,k]
        elseif  BC_xplus==BC_WALL
            u_f[Nx+1,j,k]=0.0
        end
        #水平方向
        if j!=Ny+2
            if  BC_xplus==BC_OUTLET || BC_xplus==BC_NONE
                v_f[Nx+2,j,k]=v_f[Nx+1,j,k]
            elseif  BC_xplus==BC_WALL
                v_f[Nx+2,j,k]=2*v_BC[Nx+2,j,k]-v_f[Nx+1,j,k]
            end
        end     
        if k!=Nz+2
            if  BC_xplus==BC_OUTLET || BC_xplus==BC_NONE
                w_f[Nx+2,j,k]=w_f[Nx+1,j,k]
            elseif  BC_xplus==BC_WALL
                w_f[Nx+2,j,k]=2*w_BC[Nx+2,j,k]-w_f[Nx+1,j,k]
            end
        end
    end
    #Y
    for i in 1:Nx+2,k in 1:Nz+2
        #Y-
        BC_yminus=BC_type[i,1,k]
        #垂直方向
        if BC_yminus==BC_INLET
            v_f[i,1,k]=v_BC[i,1,k]
        elseif  BC_yminus==BC_OUTLET || BC_yminus==BC_NONE
            v_f[i,1,k]=v_f[i,2,k]
        elseif  BC_yminus==BC_WALL
            v_f[i,1,k]=0.0
        end
        #水平方向
        if i!=Nx+2
            if  BC_yminus==BC_OUTLET || BC_yminus==BC_NONE
                u_f[i,1,k]=u_f[i,2,k]
            elseif  BC_yminus==BC_WALL
                u_f[i,1,k]=2*u_BC[i,1,k]-u_f[i,2,k]
            end
        end     
        if k!=Nz+2
            if  BC_yminus==BC_OUTLET || BC_yminus==BC_NONE
                w_f[i,1,k]=w_f[i,2,k]
            elseif  BC_yminus==BC_WALL
                w_f[i,1,k]=2*w_BC[i,1,k]-w_f[i,2,k]
            end
        end

        #Y+
        BC_yplus=BC_type[i,Ny+2,k]
        #垂直方向
        if BC_yplus==BC_INLET
            v_f[i,Ny+1,k]=v_BC[i,Ny+1,k]
        elseif  BC_yplus==BC_OUTLET || BC_yplus==BC_NONE
            v_f[i,Ny+1,k]=v_f[i,Ny,k]
        elseif  BC_yplus==BC_WALL
            v_f[i,Ny+1,k]=0.0
        end
        #水平方向
        if i!=Nx+2
            if  BC_yplus==BC_OUTLET || BC_yplus==BC_NONE
                u_f[i,Ny+2,k]=u_f[i,Ny+1,k]
            elseif  BC_yplus==BC_WALL
                u_f[i,Ny+2,k]=2*u_BC[i,Ny+2,k]-u_f[i,Ny+1,k]
            end
        end     
        if k!=Nz+2
            if  BC_yplus==BC_OUTLET || BC_yplus==BC_NONE
                w_f[i,Ny+2,k]=w_f[i,Ny+1,k]
            elseif  BC_yplus==BC_WALL
                w_f[i,Ny+2,k]=2*w_BC[i,Ny+2,k]-w_f[i,Ny+1,k]
            end
        end
    end
    #Z
    for i in 1:Nx+2,j in 1:Ny+2
        #Z-
        BC_zminus=BC_type[i,j,1]
        #垂直方向
        if BC_zminus==BC_INLET
            w_f[i,j,1]=w_BC[i,j,1]
        elseif  BC_zminus==BC_OUTLET || BC_zminus==BC_NONE
            w_f[i,j,1]=w_f[i,j,2]
        elseif  BC_zminus==BC_WALL
            w_f[i,j,1]=0.0
        end
        #水平方向
        if i!=Nx+2
            if  BC_zminus==BC_OUTLET|| BC_zminus==BC_NONE
                u_f[i,j,1]=u_f[i,j,2]
            elseif  BC_zminus==BC_WALL
                u_f[i,j,1]=2*u_BC[i,j,1]-u_f[i,j,2]
            end
        end     
        if j!=Ny+2
            if  BC_zminus==BC_OUTLET|| BC_zminus==BC_NONE
                v_f[i,j,1]=v_f[i,j,2]
            elseif  BC_zminus==BC_WALL
                v_f[i,j,1]=2*v_BC[i,j,1]-v_f[i,j,2]
            end
        end

        #Z+
        BC_zplus=BC_type[i,j,Nz+2]
        #垂直方向
        if BC_zplus==BC_INLET
            w_f[i,j,Nz+1]=w_BC[i,j,Nz+1]
        elseif  BC_zplus==BC_OUTLET|| BC_zplus==BC_NONE
            w_f[i,j,Nz+1]=w_f[i,j,Nz]
        elseif  BC_zplus==BC_WALL
            w_f[i,j,Nz+1]=0.0
        end
        #水平方向
        if i!=Nx+2
            if  BC_zplus==BC_OUTLET|| BC_zplus==BC_NONE
                u_f[i,j,Nz+2]=u_f[i,j,Nz+1]
            elseif  BC_zplus==BC_WALL
                u_f[i,j,Nz+2]=2*u_BC[i,j,Nz+2]-u_f[i,j,Nz+1]
            end
        end     
        if j!=Ny+2
            if  BC_zplus==BC_OUTLET|| BC_zplus==BC_NONE
                v_f[i,j,Nz+2]=v_f[i,j,Nz+1]
            elseif  BC_zplus==BC_WALL
                v_f[i,j,Nz+2]=2*v_BC[i,j,Nz+2]-v_f[i,j,Nz+1]
            end
        end
    end
end

function main()
    #パラメータ設定
    Nx=200
    Ny=200
    Nz=1

    dx=0.005
    dy=0.005
    dz=0.005

    x_min=0.0
    y_min=0.0
    z_min=0.0

    dt=0.001
    max_iter=60000

    rho=1000.0
    mu=1.0

    #変数
    pres    =zeros(Nx, Ny, Nz)
    u       =zeros(Nx+1, Ny+2, Nz+2)
    v       =zeros(Nx+2, Ny+1, Nz+2)
    w       =zeros(Nx+2, Ny+2, Nz+1)
    matA=spzeros(Nx*Ny*Nz,Nx*Ny*Nz)
    b=zeros(Nx*Ny*Nz)

    u_temp      =zeros(Nx+1, Ny+2, Nz+2)
    v_temp      =zeros(Nx+2, Ny+1, Nz+2)
    w_temp      =zeros(Nx+2, Ny+2, Nz+1)

    BC_type =fill(BC_NONE, Nx+2, Ny+2, Nz+2)
    u_BC    =zeros(Nx+1, Ny+2, Nz+2)
    v_BC    =zeros(Nx+2, Ny+1, Nz+2)
    w_BC    =zeros(Nx+2, Ny+2, Nz+1)

    #高速化用変数
    dx2_inv=1/dx/dx
    dy2_inv=1/dy/dy
    dz2_inv=1/dz/dz
    dx_inv=1/dx
    dy_inv=1/dy
    dz_inv=1/dz
    nu=mu/rho

    plts = Plots.plot()

    #初期条件設定
    for i in 1:Nx+2,j in 1:Ny+2,k in 1:Nz+2
        if i!=Nx+2
            if j==Ny+2
                u_BC[i,j,k]=1.0
            else
                u_BC[i,j,k]=0.0
            end
        end

        if j!=Ny+2
            v_BC[i,j,k]=0.0
        end

        if k!=Nz+2
            w_BC[i,j,k]=0.0
        end
    end

    #境界条件
    for i in 1:Nx+2,j in 1:Ny+2
        BC_type[i,j,1   ]=BC_NONE
        BC_type[i,j,Nz+2]=BC_NONE
    end
    for i in 1:Nx+2,k in 1:Nz+2
        BC_type[i,1   ,k]=BC_WALL
        BC_type[i,Ny+2,k]=BC_WALL
    end
    for j in 1:Ny+2,k in 1:Nz+2
        BC_type[1   ,j,k]=BC_WALL
        BC_type[Nx+2,j,k]=BC_WALL
    end

    #初期境界条件
    current_time=0
    bound_vel!(u,v,w,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min,BC_type,u_BC,v_BC,w_BC)
    #シミュレーション開始
    for num_iter in 1:max_iter
        
        current_time+=dt
        #仮速度の初期化
        u_temp=u
        v_temp=v
        w_temp=w

        #移流項・拡散項計算
        for i in 2:Nx+1,j in 2:Ny+1,k in 2:Nz+1
            #u
            if i!=Nx+1
                u_xlo1=0.5*(u[i-1,j,k]+u[i,j,k])
                u_xhi1=0.5*(u[i,j,k]+u[i+1,j,k])

                u_xlo2= 0.5*(u[i-1,j,k] + u[i,j,k])#u_xlo1>0 ? u[i-1,j,k] : u[i,j,k]
                u_xhi2= 0.5*(u[i,j,k] + u[i+1,j,k])#u_xhi1>0 ? u[i,j,k] : u[i+1,j,k]

                uu=(u_xhi1*u_xhi2-u_xlo1*u_xlo2)*dx_inv

                v_ylo= 0.5*(v[i-1,j-1,k]+v[i,j-1,k])
                v_yhi= 0.5*(v[i-1,j,k]+v[i,j,k])
                u_ylo= 0.5*(u[i,j-1,k] + u[i,j,k])#v_ylo>0 ? u[i,j-1,k] : u[i,j,k]
                u_yhi= 0.5*(u[i,j,k]   + u[i,j+1,k])#v_yhi>0 ? u[i,j,k]   : u[i,j+1,k]

                uv=(u_yhi*v_yhi-u_ylo*v_ylo)*dy_inv

                w_zlo= 0.5*(w[i-1,j,k-1]+w[i,j,k-1])
                w_zhi= 0.5*(w[i-1,j,k]+w[i,j,k])
                u_zlo= 0.5*(u[i,j,k-1] + u[i,j,k])#w_zlo>0 ? u[i,j,k-1] : u[i,j,k]
                u_zhi= 0.5*(u[i,j,k]   + u[i,j,k+1])#w_zhi>0 ? u[i,j,k]   : u[i,j,k+1]

                uw=(u_zhi*w_zhi-u_zlo*w_zlo)*dz_inv

                u_temp[i,j,k]+=-(uu+uv+uw)*dt
                #拡散項
                u_temp[i,j,k]+=nu*(u[i-1,j,k]-2*u[i,j,k]+u[i+1,j,k])*dx2_inv*dt
                u_temp[i,j,k]+=nu*(u[i,j-1,k]-2*u[i,j,k]+u[i,j+1,k])*dy2_inv*dt
                u_temp[i,j,k]+=nu*(u[i,j,k-1]-2*u[i,j,k]+u[i,j,k+1])*dz2_inv*dt
            end

            #v
            if j!=Ny+1
                u_xlo=0.5*(u[i-1,j,k]+u[i-1,j+1,k])
                u_xhi=0.5*(u[i,j,k] +u[i,j+1,k])
                v_xlo= 0.5*(v[i-1,j,k] + v[i,j,k])#u_xlo>0 ? v[i-1,j,k] : v[i,j,k]
                v_xhi= 0.5*(v[i,j,k] + v[i+1,j,k])#u_xhi>0 ? v[i,j,k] : v[i+1,j,k]

                vu=(u_xhi*v_xhi-u_xlo*v_xlo)*dx_inv

                v_ylo1= 0.5*(v[i,j,k]+v[i,j-1,k])
                v_yhi1= 0.5*(v[i,j+1,k]+v[i,j,k])
                v_ylo2= 0.5*(v[i,j-1,k] + v[i,j,k])#v_ylo1>0 ? v[i,j-1,k] : v[i,j,k]
                v_yhi2= 0.5*(v[i,j,k]   + v[i,j+1,k])#v_yhi1>0 ? v[i,j,k]   : v[i,j+1,k]

                vv=(v_yhi1*v_yhi2-v_ylo1*v_ylo2)*dy_inv

                w_zlo= 0.5*(w[i,j-1,k-1]+w[i,j,k-1])
                w_zhi= 0.5*(w[i,j-1,k]+w[i,j,k])
                v_zlo= 0.5*(v[i,j,k-1] + v[i,j,k])#w_zlo>0 ? v[i,j,k-1] : v[i,j,k]
                v_zhi= 0.5*(v[i,j,k]   + v[i,j,k+1])#w_zhi>0 ? v[i,j,k]   : v[i,j,k+1]

                vw=(v_zhi*w_zhi-v_zlo*w_zlo)*dz_inv

                v_temp[i,j,k]+=-(vu+vv+vw)*dt
                #拡散項
                v_temp[i,j,k]+=nu*(v[i-1,j,k]-2*v[i,j,k]+v[i+1,j,k])*dx2_inv*dt
                v_temp[i,j,k]+=nu*(v[i,j-1,k]-2*v[i,j,k]+v[i,j+1,k])*dy2_inv*dt
                v_temp[i,j,k]+=nu*(v[i,j,k-1]-2*v[i,j,k]+v[i,j,k+1])*dz2_inv*dt
            end

            #w
            if k!=Nz+1
                u_xlo=0.5*(u[i-1,j,k]+u[i-1,j,k+1])
                u_xhi=0.5*(u[i,j,k] +u[i,j,k+1])
                w_xlo= 0.5*(w[i-1,j,k] + w[i,j,k])#u_xlo>0 ? w[i-1,j,k] : w[i,j,k]
                w_xhi= 0.5*(w[i,j,k] + w[i+1,j,k])#u_xhi>0 ? w[i,j,k] : w[i+1,j,k]

                wu=(u_xhi*w_xhi-u_xlo*w_xlo)*dx_inv

                v_ylo= 0.5*(v[i,j-1,k-1]+v[i,j-1,k])
                v_yhi= 0.5*(v[i,j,k-1]  +v[i,j,k])
                w_ylo= 0.5*(w[i,j-1,k] + w[i,j,k])#v_ylo>0 ? w[i,j-1,k] : w[i,j,k]
                w_yhi= 0.5*(w[i,j,k]   + w[i,j+1,k])#v_yhi>0 ? w[i,j,k]   : w[i,j+1,k]

                wv=(v_yhi*w_yhi-v_ylo*w_ylo)*dy_inv

                w_zlo1= 0.5*(w[i,j,k]+w[i,j,k-1])
                w_zhi1= 0.5*(w[i,j,k+1]+w[i,j,k])
                w_zlo2= 0.5*(w[i,j,k-1] + w[i,j,k])#w_zlo>0 ? w[i,j,k-1] : w[i,j,k]
                w_zhi2= 0.5*(w[i,j,k]   + w[i,j,k+1])#w_zhi>0 ? w[i,j,k]   : w[i,j,k+1]

                vw=(w_zhi1*w_zhi2-w_zlo1*w_zlo2)*dz_inv

                w_temp[i,j,k]+=-(wu+wv+ww)*dt
                #拡散項
                w_temp[i,j,k]+=nu*(w[i-1,j,k]-2*w[i,j,k]+w[i+1,j,k])*dx2_inv*dt
                w_temp[i,j,k]+=nu*(w[i,j-1,k]-2*w[i,j,k]+w[i,j+1,k])*dy2_inv*dt
                w_temp[i,j,k]+=nu*(w[i,j,k-1]-2*w[i,j,k]+w[i,j,k+1])*dz2_inv*dt
            end
        end

        #内力項の計算(重力)u* +=(rho*g+F_s)
        #calc_internal()

        #境界条件
        bound_vel!(u_temp,v_temp,w_temp,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min,BC_type,u_BC,v_BC,w_BC)

        #圧力計算          
        for i in 1:Nx,j in 1:Ny,k in 1:Nz
            num = (i-1) + (j-1)*Nx +(k-1)*(Nx*Ny)+1
            A_diag=0.0
            if i> 1
                matA[num,num-1]=dx2_inv
                A_diag-=dx2_inv
            elseif BC_type[i,j+1,k+1]==BC_OUTLET
                A_diag-=2*dx2_inv
            end
            if i< Nx
                matA[num,num+1]=dx2_inv
                A_diag-=dx2_inv
            elseif BC_type[i+2,j+1,k+1]==BC_OUTLET
                A_diag-=2*dx2_inv
            end

            if j> 1
                matA[num,num-(Nx)]=dy2_inv
                A_diag-=dy2_inv
            elseif BC_type[i+1,j,k+1]==BC_OUTLET
                A_diag-=2*dy2_inv
            end
            if j< Ny
                matA[num,num+(Nx)]=dy2_inv
                A_diag-=dy2_inv
            elseif BC_type[i+1,j+2,k+1]==BC_OUTLET
                A_diag-=2*dy2_inv
            end

            if k< Nz
                matA[num,num+(Nx)*(Ny)]=dz2_inv
                A_diag-=dz2_inv
            elseif BC_type[i+1,j+1,k+2]==BC_OUTLET
                A_diag-=2*dyz_inv
            end
            if k> 1
                matA[num,num-(Nx)*(Ny)]=dz2_inv
                A_diag-=dz2_inv
            elseif BC_type[i+1,j+1,k]==BC_OUTLET
                A_diag-=2*dyz_inv
            end
            matA[num,num]=A_diag*1.0001
        end
        #b=zeros(Nx*Ny*Nz)
        for i in 1:Nx,j in 1:Ny,k in 1:Nz
            num = (i-1) + (j-1)*Nx +(k-1)*(Nx*Ny)+1
            divU=(u_temp[i+1,j+1,k+1]-u_temp[i,j+1,k+1])/dx
            divV=(v_temp[i+1,j+1,k+1]-v_temp[i+1,j,k+1])/dy
            divW=(w_temp[i+1,j+1,k+1]-w_temp[i+1,j+1,k])/dz
            b[num]=(divU+divV+divW)*rho/dt
        end
        answer= matA \ b
        pres=collect(reshape(answer,Nx,Ny,Nz))

        #境界条件の挿入(圧力)
        #bound_pressure()

        #速度更新
        for i in 2:Nx+1,j in 2:Ny+1,k in 2:Nz+1
            if i!=Nx+1
                u[i,j,k]=u_temp[i,j,k]-(dt/rho)*(pres[i,j-1,k-1]-pres[i-1,j-1,k-1])/dx
            end
            if j!=Ny+1
                v[i,j,k]=v_temp[i,j,k]-(dt/rho)*(pres[i-1,j,k-1]-pres[i-1,j-1,k-1])/dy
            end
            if k!=Nz+1
                w[i,j,k]=w_temp[i,j,k]-(dt/rho)*(pres[i-1,j-1,k]-pres[i-1,j-1,k-1])/dz
            end
        end

        #境界条件の挿入(流速)
        bound_vel!(u,v,w,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min,BC_type,u_BC,v_BC,w_BC)
        if num_iter%100==0
            plts=output(u,v,w,Nx,Ny,Nz,dx,dy,dz,x_min,y_min,z_min,pres)
            Plots.display(plts)
            filename = @sprintf("C:/Users/mori1/Desktop/test/%04d.png", floor(num_iter/100))
            Plots.savefig(filename)
        end
        @printf("%d\n",num_iter)
    end
end

main()

@printf("Simulation Finished")