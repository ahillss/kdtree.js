var kdtree=kdtree||(function(){
    function searchTree(P,V,invV,ray_tmax,traverseStk) {
        while(traverseStk.length>0) {
            var cur=traverseStk.pop();

            var tmin=cur.tmin;
            var tmax=cur.tmax;
            var node=cur.node;

            if(ray_tmax < tmin) {
                return null;
            }

            //leaf
            if(node.type==3) {
                return cur;
            }

            //branch
            var axis=node.type;
            var split=node.split;

            //Compute parametric distance along ray to split plane
            var tplane=(node.split-P[axis])*invV[axis];

            // Get node children pointers for ray node under plane is always
            // stored after parent node in memory aka left node is under
            var belowFirst=(P[axis]<split) || (P[axis] == split && V[axis] >= 0.0);

            //
            var firstNode=belowFirst?node.left:node.right;
            var secondNode=belowFirst?node.right:node.left;

            // Advance to next child node, possibly enqueue other child
            if(tplane > tmax || tplane <= 0.0) {
                traverseStk.push({node:firstNode,tmin:tmin,tmax:tmax});
            } else if (tplane < tmin) {
                traverseStk.push({node:secondNode,tmin:tmin,tmax:tmax});
            } else {
                traverseStk.push({node:secondNode,tmin:tplane,tmax:tmax});
                traverseStk.push({node:firstNode,tmin:tmin,tmax:tplane});
            }
        }

        return null;
    }

    //float/null <= float[3] float[3] float[3] float[3] float[3]
    function intersectAabb(P,V,invV,bMin,bMax) {
        var tnear=new Array(3);
        var tfar=new Array(3);

        for(var i=0;i<3;i++) {
            var tmin = (bMin[i]-P[i])*invV[i];
            var tmax = (bMax[i]-P[i])*invV[i];

            tnear[i] = Math.min(tmin,tmax);
            tfar[i] = Math.max(tmin,tmax);
        }

        var enter = Math.max(tnear[0], Math.max(tnear[1], tnear[2]));
        var exit = Math.min(tfar[0], Math.min(tfar[1], tfar[2]));

        return (exit > Math.max(enter, 0.0))?{tmin:enter,tmax:exit}:null;
    }

    //float <= obj float[3] float[3] float[3] float[3] float[3] func
    function intersectTree(rootNode,P,V,invV,bmin,bmax,shapeIntersect) {
		//Compute initial parametric range of ray inside kd-tree extent
        var rbRes=intersectAabb(P,V,invV,bmin,bmax);

        if(rbRes==null) {
            return Infinity;
        }

        //
        var traverseStk=[{node:rootNode,tmin:rbRes.tmin,tmax:rbRes.tmax}]

        //
        var rayMax=Infinity;
        var cur_t=Infinity;

        while(true) {
            var trv=searchTree(P,V,invV,rayMax,traverseStk);

            if(trv==null) {
                break;
            }

            //Check for intersections inside leaf node
            for(var i=0;i<trv.node.primInds.length;i++) {
                var t=shapeIntersect(P,V,invV,trv.node.primInds[i]);

                if(t!=null && t < cur_t) {
                    cur_t=t;
                }
            }

            //for shapes crossing the split plane, make sure intersection is within the current leaf's primBounds
            if(cur_t >= trv.tmin && cur_t < trv.tmax) {
                break;
            }
        }

        return cur_t;
    }

    //bool <= obj float[3] float[3] float[3] float[3] float[3] obj[] func float
    function intersectTreeP(rootNode,P,V,invV,bmin,bmax,traverseStk,shapeIntersect,rayMax) {
        var traverseNum=1;
        traverseStk[0].node=rootNode;

        var rbRes=intersectAabb(P,V,invV,bmin,bmax);

        if(rbRes==null) {
            return false;
        }

        traverseStk[0].tmin=rbRes.tmin;
        traverseStk[0].tmax=rbRes.tmax;

        //
        while(true) {
            var trv=searchTree(P,V,invV,rayMax,traverseStk,traverseNum);

            if(trv==null) {
                break;
            }

            //test against shapes within leaf
            for(var i=0;i<trv.node.primInds.length;i++) {
                var t=shapeIntersect(P,V,invV,trv.node.primInds[i]);

                if(t!=null && t <= rayMax) {
                    return true
                }
            }
        }

        return false;
    }

    //obj <= obj[] float[3] float[3] int
    function buildTree(primBounds,bmin,bmax,maxDepth) {
        var isectCost=80.0;
        var traversalCost=1.0;
        var emptyBonus=0.5;

        var maxPrims=1;
        var badRefines=0;

        maxDepth=(maxDepth!=0)?maxDepth:Math.floor(0.5+ (8.0+1.3*Math.floor(Math.log2(primBounds.length))));

        //
        var primInds=new Array(primBounds.length);

        for(var i=0;i<primBounds.length;i++) {
            primInds[i]=i;
        }

        //
        var rootNode={};

        var buildStk=[{
            node:rootNode,
            primInds:primInds,
            bmin:[bmin[0],bmin[1],bmin[2]],
            bmax:[bmax[0],bmax[1],bmax[2]],
            depth:maxDepth}]

        //
        while(buildStk.length>0) {
            var cur=buildStk.pop();

            if(cur.primInds.length <=maxPrims || cur.depth==0) {
                //Initialize leaf node if termination criteria met
                //make leaf
                cur.node.type=3;
                cur.node.primInds=cur.primInds;
            } else {
                var edges=new Array(3);

                //Choose split axis position for interior node
                var diag=[cur.bmax[0]-cur.bmin[0],cur.bmax[1]-cur.bmin[1],cur.bmax[2]-cur.bmin[2]];
                var bestAxis=-1;
                var bestOffset=-1;
                var bestCost=Infinity ;
                var oldCost=isectCost*cur.primInds.length;
                var totalSa=2.0*(diag[0]*diag[1]+diag[0]*diag[2]+diag[1]*diag[2]);//SurfaceArea
                var invTotalSa=1.0/totalSa;
                var retries=0;

                // //Choose which axis to split along
                var axis=(diag[0]>diag[1] && diag[0]>diag[2])?0:((diag[1]>diag[2])?1:2);

                // //
                do {
                    edges[axis]=new Array(cur.primInds.length*2);

                    //Initialize edges for axis
                    for(var i=0;i<cur.primInds.length;i++) {
                        var v=cur.primInds[i];
                        edges[axis][2*i  ]={type:0,val:primBounds[v].min[axis],primInd:v};
                        edges[axis][2*i+1]={type:1,val:primBounds[v].max[axis],primInd:v};
                    }

                    edges[axis].sort(function(a,b){
                        if(a.val==b.val) {
                            return (a.type < b.type)?-1:1;
                        }

                        return (a.val<b.val)?-1:1;
                    });

                    //Compute cost of all splits for axis to find best
                    var nBelow=0;
                    var nAbove=cur.primInds.length;

                    for(var i=0;i<2*cur.primInds.length;i++) {
                        if(edges[axis][i].type==1) {
                            nAbove--;
                        }

                        var edget=edges[axis][i].val;

                        //
                        if(edget>cur.bmin[axis] && edget<cur.bmax[axis]) {
                            //Compute cost for split at i th edge
                            var otherAxis0=(axis+1)%3;
                            var otherAxis1=(axis+2)%3;
                            var belowSa=2.0*(diag[otherAxis0]*diag[otherAxis1]+(edget-cur.bmin[axis])*(diag[otherAxis0]+diag[otherAxis1]));
                            var aboveSa=2.0*(diag[otherAxis0]*diag[otherAxis1]+(cur.bmax[axis]-edget)*(diag[otherAxis0]+diag[otherAxis1]));
                            var pBelow=belowSa*invTotalSa;
                            var pAbove=aboveSa*invTotalSa;
                            var eb=(nAbove==0 || nBelow==0)?emptyBonus:0.0;
                            var cost=traversalCost+isectCost*(1.0-eb)*(pBelow*nBelow+pAbove*nAbove);

                            //Update best split if this is lowest cost so far
                            if(cost < bestCost) {
                                bestCost=cost;
                                bestAxis=axis;
                                bestOffset=i;
                            }
                        }

                        if(edges[axis][i].type==0) {
                            nBelow++;
                        }
                    }

                    //
                    axis = (axis+1) % 3;

                } while(bestAxis == -1 && retries++ < 2);

                if(bestCost > oldCost) {
                    badRefines++;
                }

                if((bestCost>(4.0*oldCost) && cur.primInds.length<16) || bestAxis == -1 || badRefines == 3) {
                    //make leaf
                    cur.node.type=3;
                    cur.node.primInds=cur.primInds;
                } else {
                    //make branch
                    var split=edges[bestAxis][bestOffset].val;
                    var leftInds=[];
                    var rightInds=[];

                    // //Classify primitives with respect to split
                    for(var i=0;i<bestOffset;i++) {
                        if(edges[bestAxis][i].type==0) {
                            leftInds.push(edges[bestAxis][i].primInd);
                        }
                    }

                    for(var i=bestOffset+1;i<2*cur.primInds.length;i++) {
                        if(edges[bestAxis][i].type==1) {
                            rightInds.push(edges[bestAxis][i].primInd);
                        }
                    }

                    // //
                    cur.node.right={};
                    cur.node.left={};

                    var rightBuild={
                        node:cur.node.right,
                        primInds:rightInds,
                        bmin:[cur.bmin[0],cur.bmin[1],cur.bmin[2]],
                        bmax:[cur.bmax[0],cur.bmax[1],cur.bmax[2]],
                        depth:cur.depth-1};

                    var leftBuild={
                        node:cur.node.left,
                        primInds:leftInds,
                        bmin:[cur.bmin[0],cur.bmin[1],cur.bmin[2]],
                        bmax:[cur.bmax[0],cur.bmax[1],cur.bmax[2]],
                        depth:cur.depth-1};

                    rightBuild.bmin[bestAxis]=split;
                    leftBuild.bmax[bestAxis]=split;

                    buildStk.push(rightBuild)
                    buildStk.push(leftBuild)

                    cur.node.type=bestAxis;
                    cur.node.split=split;
                }
            }
        }

        //
        return rootNode;
    }

    return {
        "build" : buildTree,
        "search" : searchTree,
        "intersect" : intersectTree,
        "intersectP" : intersectTreeP
    };
})();
