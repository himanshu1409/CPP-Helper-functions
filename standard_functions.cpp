// Some tips
// 1. Brute Force
// 2. DP
// 3. Binary Search
// 4. Greedy
// 5. Consecutive/contiguous/subarray/substring word -> Sliding Window
// 6. Think in reverse


ll log_base2(ll n){
    ll ans=1,count=0;
    while(ans*2<=n){
        ans*=2;
        count++;
    }
    return count;
}

ll nC2(ll n){
    ll ans=(n*(n-1))/2;
    return ans;
}

ll ceil_div(ll a,ll b){
    return (a%b == 0)? a/b : a/b+1;
}

ll pwr(ll a,ll b){
    ll ans=1;
    for(int i=1;i<=b;++i)
        ans*=a;
    return ans;
}

ll pwr(ll x, ll y, ll p = MOD){
    if (y == 0)
        return 1;
    ll res = pwr(x, y / 2) % p;
    res = (res * res) % p;
    if (y % 2 == 1) res = (res * x) % p;
    return res;
}

ll pwr_without_mod(ll x, ll y){
    if (y == 0)
        return 1;
    ll res = pwr_without_mod(x, y / 2);
    res = (res * res);
    if (y % 2 == 1) res = (res * x);
    return res;
}

void primeSieve(vector<ll> &sieve,ll N=1e6){
    sieve[1] = sieve[0] = 0;
    sieve[2]=1;
    for(ll i=3; i<=N; i+=2)
        sieve[i] = 1;
    for(ll i=3; i<=N; i++){ 
        if(sieve[i]){
            for(ll j = i*i; j<=N; j = j + i)
                sieve[j] = 0;
        }
    }
}

ll power_of_2(ll n){
    ll ans=0;
    while(n%2!=0 || n==0){
        ans++;
        n/=2;
    }
    return ans;
}

ll ceil_div(ll a,ll b){
    return (a%b == 0)? a/b : a/b+1;
}

ll pwr(ll a,ll b){
    ll ans=1;
    for(int i=1;i<=b;++i)
        ans*=a;
    return ans;
}

ll factorial(ll n){
    ll prod=1;
    for(ll i=1;i<=n;++i)
        prod*=i;
    return prod;
}

bool is_power_of_2(ll n){
    return ceil(log2(n))==floor(log2(n)); // x&(x-1)==0
}

ll NcR(ll n,ll r){
    ll p = 1, k = 1;
    if (n - r < r) r = n - r;
    if (r != 0) {
        while (r) {
            p *= n;
            k *= r;
            ll m = __gcd(p, k);
            p /= m;
            k /= m;
            n--;
            r--;
        }
    }
    else
        p = 1;
    return p;
}

ll NcR_mod_m(ll n,ll r,ll m){
    ll ans=1;
    if(r==0) return n%m;
    for(ll i=1,j=n;i<=r;++i,--j)
        ans=(ans*j)%m;
    for(ll i=1;i<=r;++i)
        ans=(ans*mod_inv(i,m))%m;
    return ans;
}

ll mod_inv(ll a, ll m){
    return pwr(a, m - 2);
}

ll mod_div(ll a, ll b, ll m = M){
    return (a * mod_inv(b)) % m;
}

ll ncr(ll n, ll r, ll m=MOD){
    if(r<0 || r>n || n<0) return 0;
    return mod_div(factorial(n),(factorial(n-r)*factorial(r))%m);
}

ll maxNonOverlapIntervals(vector<pair<ll,ll>> intervals) {
    sort(intervals.begin(),intervals.end());    //Sort all intervals in ASC order
    ll count = 0;      //Count of number of intervals to be removed
    ll n = intervals.size();   //No of intervals
    ll left = 0;   //left interval
    ll right = 1;  //right interval
    
    while(right<n)  //Unless all intervals are compared
    {
        if(intervals[left].ss < intervals[right].ff)   //Non-overlapping case
        {
            left = right;
            right+=1;
        }
        else if(intervals[left].ss<=intervals[right].ss)    //Overlapping case-1 (Remove right interval)
        {
            count+=1;       //Delete right
            right+=1;
        }
        else if(intervals[left].ss>=intervals[right].ss)     //Overlapping case-2 (Remove left interval)
        {
            count+=1;
            left = right;
            right+=1;
        }
    }
    
    return n-count;
}

void permute_string(string& a,ll l,ll r,vector<string> &perm) { 
    if (l == r) perm.push_back(a);
    else{ 
        for (int i = l; i <= r; i++) { 
            swap(a[l], a[i]); 
            permute(a, l + 1, r,perm); 
            swap(a[l], a[i]); 
        } 
    } 
    return;
}

void permute_vector(vector<ll> &arr,ll l,ll r,vector<vector<ll>> &perm) { 
    if (l == r) perm.push_back(arr);
    else{ 
        for (int i = l; i <= r; i++) { 
            swap(arr[l], arr[i]); 
            permute_vector(arr, l + 1, r,perm); 
            swap(arr[l], arr[i]); 
        } 
    } 
    return;
}

ll setBits(ll n){
    return __builtin_popcount(n);
}

ll XOR_from_1_to_n(int n) { 
  if (n % 4 == 0) return n; 
  if (n % 4 == 1) return 1; 
  if (n % 4 == 2) return n + 1; 
  return 0; 
} 

ll findParent(ll node,vector<ll> &parent){
    if(parent[node]==node) return node;
    return parent[node]=findParent(parent[node],parent);
}

void union_by_rank(int u,int v,vector<int> &rank,vector<int> &parent){
    int ult_par_u=findParent(u,parent),ult_par_v=findParent(v,parent);
    if(rank[ult_par_u]>rank[ult_par_v]){
        parent[ult_par_v]=ult_par_u;
    }
    else if(rank[ult_par_v]>rank[ult_par_u]) parent[ult_par_u]=ult_par_v;
    else{
        parent[ult_par_v]=ult_par_u;
        rank[ult_par_u]++;
    }
    return;
}

class SGTree {
    vector<ll> seg;
public:
    SGTree(ll n) {
        seg.resize(4 * n + 1);
    }

    void printTree(){
        logarr(seg,0,4*n);
    }

    void build(ll ind, ll low,ll high,vector<ll> &arr) {
        if (low == high) {
            seg[ind] = arr[low];
            return;
        }

        ll mid = (low + high) / 2;
        build(2 * ind + 1, low, mid, arr);
        build(2 * ind + 2, mid + 1, high, arr);
        seg[ind] = min(seg[2 * ind + 1], seg[2 * ind + 2]);
    }

    ll query(ll ind,ll low,ll high,ll l,ll r) {
        // no overlap
        // l r low high or low high l r
        if (r < low || high < l) return INT_MAX;

        // complete overlap
        // [l low high r]
        if (low >= l && high <= r) return seg[ind];

        ll mid = (low + high) >> 1;
        ll left = query(2 * ind + 1, low, mid, l, r);
        ll right = query(2 * ind + 2, mid + 1, high, l, r);
        return min(left, right);
    }

    void update(ll ind,ll low,ll high,ll i,ll val) {
        if (low == high) {
            seg[ind] = val;
            return;
        }

        ll mid = (low + high) >> 1;
        if (i <= mid) update(2 * ind + 1, low, mid, i, val);
        else update(2 * ind + 2, mid + 1, high, i, val);
        seg[ind] = min(seg[2 * ind + 1], seg[2 * ind + 2]);
    }
};
