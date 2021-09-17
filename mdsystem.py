import numpy as np
#Simulation system of hard sphere liquid in 3D. 
#
#
#
#
#
class MDsystem(object):
    
    def __init__(self, N, L, xpos, ypos, zpos, xvel, yvel, zvel):
        #resetData()
        self.t = 0 #simulation clock
        self.N = N #num of particles
        self.L = L
        #Arrays containing positions
        self.x = xpos
        self.y = ypos
        self.z = zpos        
        #Arrays containing velocities
        self.vx = xvel
        self.vy = yvel
        self.vz = zvel
        self.time_to_hit = 1e10
        self.collision_time = np.zeros(N)
        self.partner = np.zeros(N)
        self.next_coll_particle
        self.next_partner
        #self.radius =
        self.count = 0 #number of collisions in the system
        
        for i in range(N): #set all collision times to infty at first
            self.collision_time[i] = 1e10
            
        for ip in range(N - 1): #find initial collision times for the particle pairs
            for jp in range(1, N):
                self.check_collision(ip, jp)
        
    def step(self):
        #find minimum collision time in list of collision times
        self.min_collision_time()
        #fast forward to time where collision happens
        self.move()
        self.t += self.time_to_hit
        #change velocity of the colliding particles
        self.contact()
        #reset the collision time for all particles about to collide w/ the already
        #colliding particles
        self.reset_collision_times()
        #find new collision times
        self.new_collision_times()
        self.count += 1
        
    def reset_data(self):
        self.t = 0
        
    def pbc_position(self, s, sidelength): #Position w.r.t. periodic boundary conditions
        if s > sidelength:
            s -= sidelength
        elif s < 0:
            s += sidelength
        return s
    
    def pbc_separation(self, ds, sidelength): #Distance w.r.t. periodic boundary conditions
        if ds > 0.5*sidelength:
            ds -= sidelength
        elif ds < -0.5*sidelength:
            ds += sidelength
        return ds
    
    def min_collision_time(self):
        #Sets collision time very large to find so that the min col time can be found
        for k in range(self.N):
            if self.collision_time[k] < self.time_to_hit:
                self.time_to_hit = self.collision_time[k]
                self.next_coll_particle = k
        self.next_partner = self.partner[self.next_coll_particle]
        
        
    def move(self):
        for k in range(self.N):
            self.collision_time[k] -= self.time_to_hit
            self.x[k] = self.pbc_position(self.x[k] - self.vx[k]*self.time_to_hit, self.L)
            self.y[k] = self.pbc_position(self.y[k] - self.vy[k]*self.time_to_hit, self.L)
            self.z[k] = self.pbc_position(self.z[k] - self.vz[k]*self.time_to_hit, self.L)
            
    def check_collision(self, p1, p2):
        #checks for collisions b/w particles p1 and p2 and periodic images of p2
        dvx = self.vx[p1] - self.vx[p2]
        dvy = self.vy[p1] - self.vy[p2]
        dvz = self.vz[p1] - self.vz[p2]
        v2 = dvx*dvx + dvy*dvy + dvz*dvz
        for xCell in range(-1, 2): #look at periodic images surrounding central cell
            for yCell in range(-1, 2):
                for zCell in range(-1, 2):
                    dx = self.x[p1] - self.x[p2] + xCell*self.L
                    dy = self.y[p1] - self.y[p2] + yCell*self.L
                    dz = self.z[p1] - self.z[p2] + zCell*self.L
                    bij = dx*dvx + dy*dvy + dz*dvz
                    if bij < 0:
                        r2 = dx*dx + dy*dy + dz*dz
                        disc = bij*bij - v2*(r2 - 2*self.radius)
                        if disc > 0:
                            tij = (-bij - np.sqrt(disc))/v2
                            if tij < self.collision_time[p1]:
                                self.collision_time[p1] = tij
                                self.partner[p1] = p2
                            if tij < self.collision_time[p2]:
                                self.collision_time[p2] = tij
                                self.partner[p2] = p1
        
    def contact(self): #calculates collision dynamics
        dx = self.pbc_position(self.x[self.next_coll_particle] - self.x[self.next_partner], self.L)
        dy = self.pbc_position(self.y[self.next_coll_particle] - self.y[self.next_partner], self.L)
        dz = self.pbc_position(self.z[self.next_coll_particle] - self.z[self.next_partner], self.L)
        dvx = self.vx[self.next_coll_particle] - self.vx[self.next_partner]
        dvy = self.vy[self.next_coll_particle] - self.vy[self.next_partner]
        dvz = self.vz[self.next_coll_particle] - self.vz[self.next_partner]
        b = dx*dvx + dy*dvy + dz*dvz
        delvx = -b*dx
        delvy = -b*dy
        delvz = -b*dz
        self.vx[self.next_coll_particle] += delvx
        self.vy[self.next_coll_particle] += delvy
        self.vz[self.next_coll_particle] += delvz
        self.vx[self.next_partner] -= delvx
        self.vx[self.next_partner] -= delvy
        self.vx[self.next_partner] -= delvz
        
    def reset_collision_times(self):
        #reset the coltime for all particles set to collide w/ the colliding pair
        self.collision_time[self.next_coll_particle] = 1e10
        self.collision_time[self.next_partner] = 1e10
        for k in range(self.N):
            if self.partner[k] == self.next_coll_particle:
                self.collision_time[k] = 1e10
            elif self.partner[k] == self.next_partner:
                self.collision_time[k] = 1e10
                
    def new_collision_times(self):
        #new collision time for all particles that were set to collide
        #w/ two colliding particles also finds new times for colliding particles
        for k in range(self.N):
            if k != self.next_coll_particle and k != self.next_partner:
                self.check_collision(k, self.next_partner)
                self.check_collision(k, self.next_coll_particle)