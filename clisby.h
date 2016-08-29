#include "base.h"
class Box
{
public:
	//double a.x,b.x,a.y,b.y,a.z,b.z;
	Point a,b;
	Box(){}
	Box(Point aa,Point bb){a=aa;b=bb;}
	Box(double x_min,double x_max,double y_min,double y_max,double z_min,double z_max)
	{a.x=x_min,b.x=x_max,a.y=y_min,b.y=y_max,a.z=z_min,b.z=z_max;}
	Box bound(Box bb){return Box(min(this->a.x,bb.a.x),
		max(this->b.x,bb.b.x),
		min(this->a.y,bb.a.y),
		max(this->b.y,bb.b.y),
		min(this->a.z,bb.a.z),
		max(this->b.z,bb.b.z));}
	Box intersec(Box bb){return Box(max(this->a.x,bb.a.x),
		min(this->b.x,bb.b.x),
		max(this->a.y,bb.a.y),
		min(this->b.y,bb.b.y),
		max(this->a.z,bb.a.z),
		min(this->b.z,bb.b.z));}
	Box& operator+(Point x){a+=x;b+=x;return(*this);}
	bool notempty(){if(a.x>b.x)return false;if(a.y>b.y)return false;if(a.z>b.z)return false; return true;}
	Box apply(Matrix m,Point p)
	{
		Point v2=m.dot(Point(a.x,a.y,b.z)),
			v3=m.dot(Point(a.x,b.y,a.z)),
			v4=m.dot(Point(b.x,a.y,a.z)),
			v5=m.dot(Point(b.x,b.y,a.z)),
			v6=m.dot(Point(a.x,b.y,b.z)),
			v7=m.dot(Point(b.x,a.y,b.z)),
			v1=m.dot(a),
			v8=m.dot(b);
		return Box(min(min(min(v1.x,v2.x),min(v3.x,v4.x)),min(min(v5.x,v6.x),min(v7.x,v8.x))),
			max(max(max(v1.x,v2.x),max(v3.x,v4.x)),max(max(v5.x,v6.x),max(v7.x,v8.x))),
			min(min(min(v1.y,v2.y),min(v3.y,v4.y)),min(min(v5.y,v6.y),min(v7.y,v8.y))),
			max(max(max(v1.y,v2.y),max(v3.y,v4.y)),max(max(v5.y,v6.y),max(v7.y,v8.y))),
			min(min(min(v1.z,v2.z),min(v3.z,v4.z)),min(min(v5.z,v6.z),min(v7.z,v8.z))),
			max(max(max(v1.z,v2.z),max(v3.z,v4.z)),max(max(v5.z,v6.z),max(v7.z,v8.z))))+p;
	}
	Box apply(Matrix m)
	{
		Point v2=m.dot(Point(a.x,a.y,b.z)),
			v3=m.dot(Point(a.x,b.y,a.z)),
			v4=m.dot(Point(b.x,a.y,a.z)),
			v5=m.dot(Point(b.x,b.y,a.z)),
			v6=m.dot(Point(a.x,b.y,b.z)),
			v7=m.dot(Point(b.x,a.y,b.z)),
			v1=m.dot(a),
			v8=m.dot(b);
		return Box(min(min(min(v1.x,v2.x),min(v3.x,v4.x)),min(min(v5.x,v6.x),min(v7.x,v8.x))),
			max(max(max(v1.x,v2.x),max(v3.x,v4.x)),max(max(v5.x,v6.x),max(v7.x,v8.x))),
			min(min(min(v1.y,v2.y),min(v3.y,v4.y)),min(min(v5.y,v6.y),min(v7.y,v8.y))),
			max(max(max(v1.y,v2.y),max(v3.y,v4.y)),max(max(v5.y,v6.y),max(v7.y,v8.y))),
			min(min(min(v1.z,v2.z),min(v3.z,v4.z)),min(min(v5.z,v6.z),min(v7.z,v8.z))),
			max(max(max(v1.z,v2.z),max(v3.z,v4.z)),max(max(v5.z,v6.z),max(v7.z,v8.z))));
	}
};

class tree
{
public:
	tree(){}
	tree *wl,*wr,*wp;
	double r,k;
	int n;
	Box B;
	Point Xe,X;
	double X2;
	Matrix q;
	bool full(){return (wl!=0 && wr!=0);}
};

class CClisby
{
private:
	int num2level;
	tree* root;
	tree* st;
	int nsteps; // number of steps in the walk  
	int generation;

	void Merge(tree* wl,tree* wr,tree* w)
	{
		w->r=wr->r;
		w->k=wr->k;
		w->n = wl->n+wr->n;
		w->B = wl->B.bound(wr->B.apply(w->q,wl->Xe));
		w->Xe = wl->Xe+w->q.dot(wr->Xe);
		w->X = wl->X+w->q.dot(wr->X)+wl->Xe*wr->n;
		w->X2 = wl->X2+wr->X2+wl->Xe.dot(w->q.dot(wr->X))*2+wl->Xe.dot(wl->Xe)*wr->n;
		return;
	}

	void LR(tree* w)
	{
		tree* wt;
		Matrix qt;
		wt = w->wr;
		w->wr = wt->wr;
		wt->wr = wt->wl;
		wt->wl = w->wl;
		w->wl = wt;
		qt = w->q;
		w->q = qt.dot(w->wl->q);
		w->wl->q = qt;
		Merge(w->wl->wl,w->wl->wr,w->wl);
		return;
	}
	void RR(tree* w)
	{
		tree* wt;
		Matrix qt;
		wt = w->wl;
		w->wl = wt->wl;
		wt->wl = wt->wr;
		wt->wr = w->wr;
		w->wr = wt;
		qt = w->q;
		w->q = w->wr->q;
		w->wr->q = w->q.inv().dot(qt);
		Merge(w->wr->wl,w->wr->wr,w->wr);
		return;
	}

	tree* Find_node(int nt)
	{
		int shift;
		if(nt<(root->n-num2level)) shift=num2level+nt;
		else shift=nt - (root->n-num2level);

		return st+shift;
	}

	tree* line_initialize()
	{
		generation = 0;
		int n=nsteps+1, n_node = n, tolnum=n_node+n_node-1;
		int child_count=0,parent_count=n_node;
		st = new tree[tolnum];
		int level=1;
		while (n>2)
		{
			if(n%2==1){num2level+=twopowof(level);level++;n=((n+1)/2)-1;}
			else {level++;n=n/2;}
		}

		for (int i=0;i<tolnum-1;i++)
		{
			st[i].q=I;
			if(i<n_node)
			{
				st[i].r=RADIUS;
				st[i].k=RIGID;
				st[i].wl=0;
				st[i].wr=0;
				st[i].n=1;
				if(i<num2level)
				{
					if(i%2==1)
					{
						st[i].B=Box(1-st[i].r,1+st[i].r,-st[i].r,st[i].r,-st[i].r,st[i].r);
						st[i].Xe = Point(1,0,0);
						st[i].X = Point(1,0,0);
					}
					else
					{
						st[i].B=Box(-st[i].r,st[i].r,1-st[i].r,1+st[i].r,-st[i].r,st[i].r);
						st[i].Xe = Point(0,1,0);
						st[i].X = Point(0,1,0);
					}
				}
				else
				{
					if(i%2==0)
					{
						st[i].B=Box(1-st[i].r,1+st[i].r,-st[i].r,st[i].r,-st[i].r,st[i].r);
						st[i].Xe = Point(1,0,0);
						st[i].X = Point(1,0,0);
					}
					else
					{
						st[i].B=Box(-st[i].r,st[i].r,1-st[i].r,1+st[i].r,-st[i].r,st[i].r);
						st[i].Xe = Point(0,1,0);
						st[i].X = Point(0,1,0);
					}
				}
				st[i].X2 = 1;
				if(i==num2level)
				{
					st[i].B=Box(-st[i].r,st[i].r,-st[i].r,st[i].r,-st[i].r,st[i].r);
					st[i].Xe = Point(0,0,0);
					st[i].X = Point(0,0,0);
					st[i].X2 = 0;
				}
			}
			else
			{
				st[i].wl=st+child_count;
				child_count++;
				st[i].wr=st+child_count;
				child_count++;
				Merge(st[i].wl,st[i].wr,st+i);
			}
			st[i].wp=st+parent_count;
			if(i%2==1)parent_count++;
		}

		st[tolnum-1].q=I;
		st[tolnum-1].wl=st+child_count;
		child_count++;
		st[tolnum-1].wr=st+child_count;
		child_count++;
		Merge(st[tolnum-1].wl,st[tolnum-1].wr,st+tolnum-1);
		st[tolnum-1].wp=0;

		return &st[tolnum-1];
	}

	//start from 1
	void Shuffle_up(int n0,tree* w)
	{
		if (n0==w->wl->n) return;
		else if(n0<w->wl->n)
		{
			Shuffle_up(n0,w->wl);
			RR(w);
		}
		else if(n0>w->wl->n)
		{
			Shuffle_up(n0 - w->wl->n,w->wr);
			LR(w);
		}
		return;
	}
	void Shuffle_down(tree* w)
	{
		int nt = int((w->n+1)/2);
		if(nt == w->wl->n) return;
		else if (nt<w->wl->n)
		{
			RR(w);
			Shuffle_down(w->wr);
		}
		else if (nt>w->wl->n)
		{
			LR(w);
			Shuffle_down(w->wl);
		}
		return;
	}
	bool Intersect(Point xla,Matrix qla,tree* wl,Point xra,Matrix qra,tree* wr )
	{
		Box blt = wl->B.apply(qla,xla);
		Box brt = wr->B.apply(qra,xra);
		if(blt.intersec(brt).notempty()==false) return false;
		if(wl->n<=1 && wr->n<=1) 
		{
			if(((qla.dot(wl->Xe)+xla)-(qra.dot(wr->Xe)+xra)).norm()<(wl->r+wr->r)) return true; 
			else return false;
		}
		if(wl->n>=wr->n)
		{
			if(Intersect(xla+qla.dot(wl->wl->Xe),qla.dot(wl->q),wl->wr,xra,qra,wr)) return true;
			return Intersect(xla,qla,wl->wl,xra,qra,wr);
		}
		else
		{
			if(Intersect(xla,qla,wl,xra,qra,wr->wl)) return true;
			return Intersect(xla,qla,wl,xra+qra.dot(wr->wl->Xe),qra.dot(wr->q),wr->wr);
		}
	}

	bool Shuffle_intersect(tree* w,Matrix q0,int wlc, int ilc)
	{
		if(wlc==1){if(Intersect(origin,I,w->wl,w->wl->Xe+w->q.dot(q0.dot(w->wr->wl->Xe)),w->q.dot(q0.dot(w->wr->q)),w->wr->wr))return true;}
		else if(wlc==0){if(Intersect(origin,I,w->wl->wl,w->wl->Xe,w->q.dot(q0),w->wr))return true;}
		else if(wlc==-1){if(Intersect(origin,I,w->wl,w->wl->Xe,w->q.dot(q0),w->wr))return true;}

		if(w->wp == 0) return false;
		int ilc_new;
		if(w->wp->wp==0){ilc_new=-1;}
		else if(w->wp->wp->wl==w->wp) ilc_new=true;
		else ilc_new = false;

		tree *wt = w->wp;
		if(ilc) RR(wt);
		else LR(wt);

		return Shuffle_intersect(wt,q0,ilc,ilc_new);
	}

// 	void Pseudo_dimerize(tree* w)
// 	{
// 		if(w->wl->full())Pseudo_dimerize(w->wl);
// 		if(w->wr->full())Pseudo_dimerize(w->wr);
// 
// 		int nt;
// 		Matrix q,qt;
// 		int k;
// 		for(k=0;k<1000;k++)
// 		{
// 			if (w->wl->full())
// 			{
// 				nt = w->wl->n - Random_integer_log(1,w->wl->n);
// 				qt = Random_symmetry();
// 				Attempt_pivot(qt,nt);
// 			}
// 			if (w->wr->full())
// 			{
// 				nt = Random_integer_log(1,w->wr->n);
// 				qt = Random_symmetry();
// 				Attempt_pivot(qt,nt);
// 			}
// 			q = Random_symmetry();
// 			if(!Intersect(origin,I,w->wl,w->wl->Xe,w->q,w->wr)) break;
// 		}
// 		if(k==1000)printf("error(Pseudo_dimerize): 1000 trial not success\n");
// 
// 		Merge(w->wl,w->wr,w);
// 
// 		for (int i=1;i<=sqrt(double(w->n));i++)
// 		{
// 			nt = w->wl->n - Random_integer_log(1,w->wl->n);
// 			qt = Random_symmetry();
// 			Attempt_pivot(qt,nt);
// 			nt = w->wl->n + Random_integer_log(0,w->wr->n);
// 			qt = Random_symmetry();
// 			Attempt_pivot(qt,nt);
// 		}
// 		return;
// 	}

	tree* scan(FILE *fptr)
	{
		fscanf_s(fptr, "N=%d ind=%d ", &nsteps, &generation);
		printf("Start SAW is: N=%d ind=%d \n", nsteps, generation);
		Sphere curstep;
		Point laststep;
		laststep.zero();

		int n=nsteps+1, n_node = n, tolnum=n_node+n_node-1;
		int child_count=0,parent_count=n_node;
		st = new tree[tolnum];
		int level=1;
		while (n>2)
		{
			if(n%2==1){num2level+=twopowof(level);level++;n=((n+1)/2)-1;}
			else {level++;n=n/2;}
		}
		
		tree* mostleft = st+num2level;
		for (int i=0;i<n_node;i++)
		{
			int j = (i<(n_node-num2level)?i:(i-n_node));
			(mostleft+j)->wl=0;
			(mostleft+j)->wr=0;
			(mostleft+j)->n=1;

			curstep.scan(fptr);
			(mostleft+j)->r=curstep.r;
			(mostleft+j)->k=curstep.k;

			(mostleft+j)->Xe = curstep.center()-laststep;
			(mostleft+j)->X = (mostleft+j)->Xe;
			(mostleft+j)->X2 = (mostleft+j)->Xe.dot((mostleft+j)->Xe);
			if((mostleft+j)->Xe.dot((mostleft+j)->Xe)!=1)printf("%lf \n",(mostleft+j)->Xe.dot((mostleft+j)->Xe));
			(mostleft+j)->B=Box((mostleft+j)->Xe.x-(mostleft+j)->r, (mostleft+j)->Xe.x+(mostleft+j)->r, (mostleft+j)->Xe.y-(mostleft+j)->r, (mostleft+j)->Xe.y+(mostleft+j)->r, (mostleft+j)->Xe.z-(mostleft+j)->r, (mostleft+j)->Xe.z+(mostleft+j)->r);
			laststep = curstep.center();
		}
		for (int i=0;i<tolnum-1;i++)
		{
			st[i].q=I;
			if(i<n_node)
			{
			}
			else
			{
				st[i].wl=st+child_count;
				child_count++;
				st[i].wr=st+child_count;
				child_count++;
				Merge(st[i].wl,st[i].wr,st+i);
			}
			st[i].wp=st+parent_count;
			if(i%2==1)parent_count++;
		}
		st[tolnum-1].q=I;
		st[tolnum-1].wl=st+child_count;
		child_count++;
		st[tolnum-1].wr=st+child_count;
		child_count++;
		Merge(st[tolnum-1].wl,st[tolnum-1].wr,st+tolnum-1);
		st[tolnum-1].wp=0;

		return &st[tolnum-1];
	}

	double GetRg2()
	{
		Point rc = GetStepi(0).center();
		double dis = 0;
		double rg2 = 0;

		for (int i=1;i<=nsteps;i++)
		{
			rc = rc + GetStepi(i).center();
		} 
		rc /= (double)(nsteps+1);

		for (int i=0;i<=nsteps;i++)
		{
			dis = (GetStepi(i).center()-rc).norm();
			rg2 += dis*dis;
		} 
		rg2 /= (nsteps+1);
		return(rg2);
	}
	void deallocate()
	{
		delete[] st;
	}
public:
	CClisby(){}

	// if name is "0", initialize SAW with step n, if not, read n from file
	void ImportSAW(const char* name,int n=100)
	{
		num2level=0;
		nsteps = n;
		//if (name[0] == '0') line_initialize(2);
		if (name == "0") root = line_initialize(); //Pseudo_dimerize(w);
		else
		{
			FILE *fptr = fopen(name, "r");
			if (fptr == NULL)
			{
				printf("FAILURE to open walk file %s, exiting \n", name);
				exit(0);
			}
			root = this->scan(fptr);
			fclose(fptr);
		}
	}

	bool Attempt_pivot(Matrix qt,int nt)
	//bool Attempt_pivot_simple(Matrix qt,int nt)
	{
		Shuffle_up(nt+1,root); //start from 1
		root->q = root->q.dot(qt);
		bool intersection=Intersect(origin,I,root->wl,root->wl->Xe,root->q,root->wr);
		if(intersection) root->q = root->q.dot(qt.inv());
		else Merge(root->wl,root->wr,root);
		Shuffle_down(root);
		return (!intersection);
	}

	//bool Attempt_pivot(Matrix qt,int nt)
	bool Attempt_pivot_fast(Matrix qt,int nt) // problem in shuffle down
	{
		tree* wt = Find_node(nt);
		wt = wt->wp;
		bool ilc;
		if(wt->wp->wl==0) {ilc=true;}
		else if(wt->wp->wl==wt){ilc=true;}
		else ilc=false;
		
		bool intersection = Shuffle_intersect(wt,qt,-1,ilc);
		Shuffle_down(root);
		if(!intersection)
		{
			Shuffle_up(nt+1,root); //start from 1
			root->q = root->q.dot(qt);
			Shuffle_down(root);
			Merge(root->wl,root->wr,root);
		}
		return (!intersection);
	}

	Sphere GetStepi(int i)
	{
		tree *wlt, *wt = Find_node(i);
		double r = wt->r, k = wt->k;
		Point pt = wt->Xe;
		
		while(wt->wp != 0)
		{
			wlt = wt->wp->wl;
			if (wlt!=wt) pt = wt->wp->q.dot(pt)+wlt->Xe;
			wt = wt->wp;
		}
		return Sphere(pt.x,pt.y,pt.z,r,k);
	}
	int GetNSteps(){return nsteps;}
	void IncreaseGeneration(){generation++;}

	void Writedown(const char* name)
	{
		FILE *fptr;
		char buffer[50]; // <- danger, only storage for 256 characters.
		sprintf_s(buffer, "%s%d", name,generation);

		fptr = fopen(buffer, "w");
		fprintf(fptr, "N=%d ind=%d \n", nsteps, generation);
		for (int i = 0; i <= nsteps; i++) GetStepi(i).print(fptr); //one step less
		fclose(fptr);
	}
	void Record(const char* name)
	{
		FILE *fptr;
		char buffer[50]; // <- danger, only storage for 256 characters.
		sprintf_s(buffer, "%s_%d", name,nsteps);
		// record the "data"
		double endnorm = root->Xe.norm();
		fptr = fopen(buffer, "a");
		fprintf(fptr,"%14.10f    %14.10f \n",GetRg2(),endnorm*endnorm);
		fclose(fptr);
	}
};