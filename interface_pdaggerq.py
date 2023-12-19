import pdaggerq
import numpy as np
import re
import sys
# Replace by the correct path
sys.path.append('#PATH_TO_SimpleWick')
import copy
import wick as w
from IPython.display import display, Latex

# Factorial of a number
def factorial(k):
    if type(k) != type(1):
        print('Error type arg k in factorial')
        sys.exit()
        
    f = 1
    for i in range(2,k+1):
        f = f * i
        
    return f

def gen_left_str(ops):
    if type(ops) != type(['aa','bb']):
        print('Error type arg ops in gen_left_str')
        sys.exit()
        
    acc = 'e'+str(len(ops)//2)+'('
    for op in ops:
        acc = acc + op + ',' 
    acc = acc[0:len(acc)-1] + ')'
    acc = acc.replace('-','')
    acc = acc.replace('+','')
    return acc

# Class for a single T operator
class T_pq():
    def __init__(self,t):
        # Init
        self.is_ok = True
        self.t = t
        self.na_i = 0
        self.nb_i = 0
        self.ng_i = 0
        self.na_a = 0
        self.nb_a = 0
        self.ng_a = 0
        self.sign = 1
    
    # clean the T
    def clean(self):
        T = self.t
        #print('Not cleaned T:',T)
        T = T.replace('t','')
        T = re.sub(r'[()]','',T)
        T = re.sub(r'[0-9]+', '', T)
        T = T.split(',')
        self.t = T
        #print('Cleaned T:',T)

    # t2 act -> act = 0
    # t1 act -> act = 0
    def remove_act_exc(self):
        acc = []
        for op in self.t:
            acc.append(op[1])
        nb_n = acc.count('n')
        nb_m = acc.count('m')
        # Remove the t1 act -> act
        if len(self.t) == 2:
            if nb_m == 1 and nb_n == 1:
                self.is_ok = True
        # Remove the t2 act -> act
        elif len(self.t) == 4:
            if nb_m == 2 and nb_n == 2:
                self.is_ok = False
        
    # Count alpha, beta and general spin
    def count_spin(self):
        na_i = self.na_i
        nb_i = self.nb_i
        ng_i = self.ng_i
        na_a = self.na_a
        nb_a = self.nb_a
        ng_a = self.ng_a
        
        for op in self.t:
            if (op[2] == 'a' and op[0] == 'i'): 
                na_i = na_i + 1
            if (op[2] == 'b' and op[0] == 'i'): 
                nb_i = nb_i + 1
            if (op[2] == 'g' and op[0] == 'i'): 
                ng_i = ng_i + 1
            if (op[2] == 'a' and op[0] == 'a'): 
                na_a = na_a + 1
            if (op[2] == 'b' and op[0] == 'a'): 
                nb_a = nb_a + 1
            if (op[2] == 'g' and op[0] == 'a'): 
                ng_a = ng_a + 1
        self.na_i = na_i
        self.nb_i = nb_i
        self.ng_i = ng_i
        self.na_a = na_a
        self.nb_a = nb_a
        self.ng_a = ng_a
        #print('Spin i:',self.na_i,self.nb_i,self.ng_i)
        #print('Spin a:',self.na_a,self.nb_a,self.ng_a)
        
        # Check is the term is zero by spin
        # alpha -> + 1, beta -> -1
        s_i = na_i - nb_i
        s_a = na_a - nb_a
        if (s_i != s_a) and (abs(s_i - s_a) != ng_i + ng_a):
            self.is_ok = False
            
        #print('is ok:', self.is_ok)
        
    # Put the hole indexes before the particle ones
    def to_std_order(self):
        acc = []
        for i in range(len(self.t)//2,len(self.t)):
            acc.append(self.t[i])
        for i in range(0,len(self.t)//2):
            acc.append(self.t[i])
        #print('acc',acc)
        self.t = acc

    # Remove the first index that define the orbital class
    def remove_first_idx(self):
        acc = []
        for t in self.t:
            acc.append(t[1:3])
        self.t = acc       
        
    # Move indexes to end up with: beta g alpha
    def move_b_to_left(self):
        if len(self.t) == 2:
            return
        if len(self.t) > 4:
            print('Error, only done for t1 and t2')
            sys.exit()
        sign = 1
        idx_spin = 2
        t = copy.deepcopy(self.t)
        #print(t)
        t_i = t[0:(len(t)//2)]
        t_a = t[(len(t)//2):len(t)]
        #print("ti",t_i)
        #print("ta",t_a)
        
        # b=0, g=1, a=2
        # For hole part
        acc_i = []
        for elem in t_i:
            acc_i.append(elem[idx_spin])
        for i in range(len(acc_i)):
            if acc_i[i] == 'b':
                acc_i[i] = 0
            elif acc_i[i] == 'g':
                acc_i[i] = 1
            else:
                acc_i[i] = 2
        #print(acc_i)
        if acc_i[1] < acc_i[0]:
            sign = -sign
            t[0] = t_i[1]
            t[1] = t_i[0]

        # For particle part
        acc_a = []
        for elem in t_a:
            acc_a.append(elem[idx_spin])
        for i in range(len(acc_a)):
            if acc_a[i] == 'b':
                acc_a[i] = 0
            elif acc_a[i] == 'g':
                acc_a[i] = 1
            else:
                acc_a[i] = 2
        #print(acc_a)
        if acc_a[1] < acc_a[0]:
            sign = -sign
            t[2] = t_a[1]
            t[3] = t_a[0]

        #print(self.t,t)
        #print(self.sign,sign)
        # New t and update the sign
        self.t = t
        self.sign = sign

# Class for the whole term comming from pdaggerq (prefactor + Ts)
class Term_pq():
    def __init__(self,a_term):
        self.prefactor = float(a_term[0]) # prefactor 
        self.str_T = a_term[1:] # output from pdaggerq
        self.nb_T = len(self.str_T) # Number of T
        self.l_T = None # List of T
        self.tex = None # Tex
        #self.ref = None
        #self.deltas = deltas

        # Loop over the Ts
        acc = []
        for i in range(4):
          if i >= self.nb_T:
              break
          t = T_pq(self.str_T[i])
          # Cleaning
          t.clean()
          # Remove act -> act amplitudes
          t.remove_act_exc()
          # t(i,j,...a,b,...)
          t.to_std_order()
          #print("T n°",i,":",t.t)
          # To check if the T is zero
          t.count_spin()
          # Put the beta indexes on the left and change the sign
          t.move_b_to_left()
          #print(self.prefactor,t.sign)
          self.prefactor = self.prefactor * t.sign
          #print(self.prefactor)
          # The first index for the orbital class is not useful anymore
          t.remove_first_idx()
          # Nullify if a term is zero by spin
          if not t.is_ok:
              self.prefactor = 0.0
          acc.append(t.t)

        # List of Ts
        self.l_T = acc
        #if self.prefactor != 0.0:
        #    print('Term:',self.prefactor,self.l_T)

    #def remove_disconnected(self):

    # Convert the term to latex
    def to_latex(self):
        # Sign + prefactor
        if str(self.prefactor)[0] == '-':
            acc = str(self.prefactor)
        else:
            acc = '+' + str(self.prefactor)

        # Ts
        self.tex = acc + Ts_to_tex(self.l_T)
        #for t in self.l_T:
        #    acc = acc + ' \\ t_{'
        #    #print('t',t)
        #    # Lower indexes
        #    for i in range(0,len(t)//2):
        #        acc = acc + t[i][0] + '$' + t[i][1]
        #    acc = acc + '}^{'
        #    # Upper indexes
        #    for i in range(len(t)//2,len(t)):
        #        acc = acc + t[i][0] + '$' + t[i][1]
        #    acc = acc + '}'
        ##print(acc)
        ## Spin
        #acc = acc.replace('$a','_{\\alpha}')
        #acc = acc.replace('$b','_{\\beta}')
        #acc = acc.replace('$g','_{}')
        #self.tex = acc

    # To diplay the latex eq
    def tex_show(self):
        null = self.to_latex()
        display(Latex(f'${self.tex}$'))

def prefactor_to_tex(prefactor):
    if type(prefactor) != type(1.0):
        print('Error type arg prefactor in prefactor_to_tex')
        sys.exit()
    # Sign + prefactor
    if str(prefactor)[0] == '-':
        tex = str(prefactor)
    else:
        tex = '+' + str(prefactor)
    return tex

def Ts_to_tex(Ts):
    if type(Ts) != type([['aa'],['bb']]):
        print('Error type Ts in Ts_to_tex')
        sys.exit()
        
    tex = ''
    for t in Ts:
        tex = tex + ' \\ t_{'
        #print('t',t)
        # Lower indexes
        for i in range(0,len(t)//2):
            if t[i][1] == 'b':
                tex = tex + '\\bar{' + t[i][0] +'}'
            else:
                tex = tex + t[i][0]
            #tex = tex + t[i][0] + '$' + t[i][1]
        tex = tex + '}^{'
        # Upper indexes
        for i in range(len(t)//2,len(t)):
            if t[i][1] == 'b':
                tex = tex + '\\bar{' + t[i][0] +'}'
            else:
                tex = tex + t[i][0]
            #tex = tex + t[i][0] + '$' + t[i][1]
        tex = tex + '}'
    #print(tex)
    # Spin
    #tex = tex.replace('$a','_{\\alpha}')
    #tex = tex.replace('$b','_{\\beta}')
    #tex = tex.replace('$g','_{}')

    return tex

def deltas_to_tex(deltas):
    if type(deltas) != type([['aa','bb'],['aa','bb']]):
        print('Error type arg deltas in deltas to tex')
        sys.exit()
    #tex = w.deltas_to_tex(deltas)
    #tex = tex.replace('_{\\alpha}','')
    #tex = tex.replace('_{\\beta}','')
    tex = ''
    for delta in deltas:
        if delta[0][3] == 'b':
            d1 = '\\bar{'+str(delta[0][1])+'}'
        else:
            d1 = str(delta[0][1])
        if delta[1][3] == 'b':
            d2 = '\\bar{'+str(delta[1][1])+'}'
        else:
            d2 = str(delta[1][1])
        tex = tex + '\delta('+d1+','+d2+') \ '

    return tex

class T():
    def __init__(self,t,ref,list_act_idx):
        self.t = t
        self.kind = len(self.t)//2
        self.list_act_idx = list_act_idx
        self.ref = ref
        self.is_disconnected = self.check_connection()

    def check_connection(self):
        n = 0
        for idx in self.list_act_idx:
            for label in self.t:
                n = n + label.count(idx)

        if n == 0:
            res = True
        else:
            res = False
            
        return res

    def apply_permutation_t(self,list_perm):
        for i in range(len(self.t)):
            #for perm in list_perm:
            for perm in list_perm:
                label1 = perm[0]
                label2 = perm[1]
                #print(self.t[i],label1,label2,spin1,spin2)
                if self.t[i][0] == label1:
                    self.t[i] = label2 + self.t[i][1]
                elif self.t[i][0] == label2:
                    self.t[i] = label1 + self.t[i][1]
                #print('a',self.t[i])

class Term():
    def __init__(self,deltas,prefactor,Ts):
        #print('d:',deltas)
        self.deltas = deltas
        self.prefactor = prefactor
        self.Ts = Ts
        self.nb_T = len(Ts)
        self.is_disconnected = False
        for t in self.Ts:
            if t.is_disconnected:
                self.is_disconnected = True

    def apply_permutation_term(self,sign,list_perm):
        # Prefactor
        self.prefactor = self.prefactor * sign
        # Ts
        for t in self.Ts:
            #print("t b",t.t)
            t.apply_permutation_t(list_perm)
            
        # Delta
        for perm in list_perm:
            label1 = perm[0]
            label2 = perm[1]
            #print(perm,spin)
            #print('d b',self.deltas)
            for i in range(len(self.deltas)):
                for j in range(2):
                    if self.deltas[i][j][0] == label1:
                        self.deltas[i][j] = label2 + self.deltas[i][j][1]
                    elif self.deltas[i][j][0] == label2:
                        self.deltas[i][j] = label1 + self.deltas[i][j][1]           
            #print('d a',self.deltas)
                
    def spin_flip(self,ref_to_flip,res_ref):
        for i in range(len(self.Ts)):
            if self.Ts[i].ref != ref_to_flip:
                continue
            else:
                self.Ts[i].ref = res_ref
                
            for j in range(len(self.Ts[i].t)):
                if self.Ts[i].t[j][1] == 'a':
                    self.Ts[i].t[j] = self.Ts[i].t[j][0]+'b'
                elif self.Ts[i].t[j][1] == 'b':
                    self.Ts[i].t[j] = self.Ts[i].t[j][0]+'a'
                else:
                    print('Unknow spin for spin flip')
                    sys.exit()
        
    def Ts_to_fortran(self,shift):
        code = ''
        for t in self.Ts:
            if len(code) > 60:
                code = code + ' & \n' + shift
            tmp = str(t.t).replace('\'','')
            tmp = tmp.replace('[','(')
            tmp = tmp.replace(']',')')
            tmp = 't' + str(t.kind) + '_' + t.ref + tmp
            code = code + ' * ' + tmp
        code  = code + ' & \n'
    
        return code
    
    def deltas_to_tex(self):
        tex = deltas_2_tex(self.deltas)

        return tex

    def prefactor_to_tex(self):
        # Prefactor
        if self.prefactor == 1.0:
            tex = '+'
        elif self.prefactor == -1.0:
            tex = '-'
        elif self.prefactor > 1.0:
            tex == '+' + str(self.prefactor)
        elif self.prefactor < -1.0:
            tex == '-' + str(self.prefactor)
        else:
            print("What ? ", str(self.prefactor))
            sys.exit()

        return tex

    def Ts_to_tex(self):
        tex = ''
        # Ts
        for t in self.Ts:
            tex = tex + '\\ ^{'+t.ref+'}t_{'
            
            ## Lower indexes
            for i in range(0,len(t.t)//2):
                if t.t[i][1] == 'b':
                    tex = tex + '\\bar{' + t.t[i][0] +'}'
                else:
                    tex = tex + t.t[i][0]
                 
            ## Upper indexes
            tex = tex + '}^{'
            for i in range(len(t.t)//2,len(t.t)):
                if t.t[i][1] == 'b':
                    tex = tex + '\\bar{' + t.t[i][0] +'}'
                else:
                    tex = tex + t.t[i][0]
            tex = tex + '}'

        return tex

    def to_tex(self):
        tex = self.prefactor_to_tex()
        tex = tex + self.deltas_to_tex()
        tex = tex + self.Ts_to_tex()

        return tex

    def to_tex_no_delta(self):
        tex = self.prefactor_to_tex()
        tex = tex + self.Ts_to_tex()
        return tex

    def eq_show(self):
        txt = self.prefactor + ' ' +self.deltas + ' ' + self.Ts.t

    def tex_show(self):
        tex = self.to_tex()
        null = display_tex(tex)
        
def delta_2_tex(delta):
    tex = ''
    if delta[0][1] == 'b':
        d1 = '\\bar{'+delta[0][0]+'}'
    else:
        d1 = delta[0][0]
    if delta[1][1] == 'b':
        d2 = '\\bar{'+delta[1][0]+'}'
    else:
        d2 = delta[1][0]
        
    tex = tex + '\\delta(' + d1 + ',' + d2 + ') \ '

    return tex

def deltas_2_tex(deltas):
    tex = ''
    for delta in deltas:
        tex = tex + delta_2_tex(delta)

    return tex

class LTerms():
    def __init__(self):
        self.terms = []
        self.unique_deltas = []
        self.factorized_terms = []
        self.factorized = False

    def append_Term(self,term1):
        self.terms.append(term1)

    def append_LTerms(self,lterms1):
        for term1 in lterms1.terms:
            self.terms.append(term1)

    def append_prod(self,sign,lterms1,lterms2):
        prod = prod_LTerms(self,sign,lterms1,lterms2)
        self.append_LTerms(prod)

    def prod_LTerms(self,sign,lterms1,lterms2):
        for term1 in lterms1.terms:
            #print("1",str(term1.deltas),term1.Ts_to_fortran("")[:-4])
            for term2 in lterms2.terms:
                #print("2    ",str(term2.deltas),term2.Ts_to_fortran("")[:-4])
                
                # Product of the kronecker deltas
                is_ok = True
                if len(term1.deltas) != 0:
                    deltas = copy.deepcopy(term1.deltas)
                    for d2 in term2.deltas:
                        for d in term1.deltas:
                            #print(d,d2)
                            # check is there is two time the same operators in the deltas
                            is_ok = not(is_conflict_deltas(d,d2))
                            if not(is_ok):
                                break
                        if not(is_ok):
                            break
                        else:
                            #print("append",d2)
                            deltas.append(d2)
                else:
                    deltas = copy.deepcopy(term2.deltas)
                if not(is_ok):
                    continue

                #print("add",deltas)
                # Product of prefactors
                prefactor = sign * term1.prefactor * term2.prefactor

                # Product of Ts
                Ts = copy.deepcopy(term1.Ts)
                for t2 in term2.Ts:
                    Ts.append(t2)

                self.append_Term(Term(deltas,prefactor,Ts))

    def factorize(self):
        self.extract_unique_deltas()
        self.factorized_terms = [[] for i in range(len(self.unique_deltas))]
        for term in self.terms:
            idx = search_idx(term.deltas,self.unique_deltas)
            self.factorized_terms[idx].append(Term([],term.prefactor,term.Ts))

    def show_tex_factorized(self):
        k = 0
        for deltas in self.unique_deltas:
            tex = ""
            #if len(deltas) == 0:
            #    tex = "&"
            #else:
            #    tex = tex + "\\\ &+ "
            tex = tex + deltas_2_tex(deltas)
            tex = tex + '\\bigl['
            l = 0
            for term in self.factorized_terms[k]:
                tex = tex + term.to_tex_no_delta()
                #l += len(term.to_tex_no_delta())
                #if l > 300:
                #    tex = tex + "\\\ & "
                #    l = 0
            #if l == 0:
            #    tex = tex[:-5]
            tex = tex + '\\bigr]'
            null = display_tex(tex)
            #print(tex)
            k = k + 1
            
    def extract_unique_deltas(self):
        unique_deltas = []
        for term in self.terms:
            if not(is_in(term.deltas,unique_deltas)):
                unique_deltas.append(term.deltas)

        # Sort depending on the number of deltas
        ## Max number
        max_len = 0
        for deltas in unique_deltas:
            if len(deltas) > max_len:
                max_len = len(deltas)

        ## Split depending on the length
        acc = [[] for i in range(max_len+1)]
        for deltas in unique_deltas:
            acc[len(deltas)].append(deltas)

        ## Reduction of the number of dimensions
        tmp = []
        for list_elem in acc:
            for elem in list_elem:
                tmp.append(elem)

        self.unique_deltas = tmp

    # Bad function for bad things...
    def reverse_deltas_order(self):
        for term in self.terms:
            if len(term.deltas) >= 2:
                stop = 0
                while term.deltas[stop][1][0] == 'a' or  term.deltas[stop][1][0] == 'b':
                    stop += 1
                    if stop == len(term.deltas)-1:
                        break
                if stop >= 1:
                    term.deltas = self.move_elements(term.deltas,stop)
    
    def move_elements(self, lst, index):
        if index < 0 or index >= len(lst):
            raise ValueError("Index out of range")

        # Move elements before the index to the end
        result = lst[index:] + lst[:index]

        return result     

    def ordered_by_t1(self):
        if len(self.factorized_terms) == 0:
            return

        acc = [[] for i in range(len(self.factorized_terms))]
        for i in range(len(self.factorized_terms)):
            for term in self.factorized_terms[i]:
                # Not a T1
                if len(term.Ts[0].t) != 2:
                    continue
                # Inactive labels
                label_h = term.Ts[0].t[0]
                label_p = term.Ts[0].t[0]
                h = (label_h == 'i' or label_h == 'j')
                p = (label_p == 'a' or label_p == 'b')
                if h and p:
                    acc[i].append(term)

            # All the remaining terms
            for term in self.factorized_terms[i]:
                if not(is_in(term,acc[i])):
                    acc[i].append(term)
                      
        self.factorized_terms = acc

    def remove_disconnected(self):
        k = 0
        for i in range(len(self.terms)):
            #null = display_tex(self.terms[k].Ts_to_tex())
            #print(self.terms[k].is_disconnected)
            if self.terms[k].is_disconnected:
                self.terms.pop(k)
            else:
                k = k + 1

    def spin_flip(self,ref_to_flip,res_ref):
        for term in self.terms:
            term.spin_flip(ref_to_flip,res_ref)

    def apply_permutation(self,sign,list_perm):
        for term in self.terms:
            term.apply_permutation_term(sign,list_perm)

        #acc = []
        #for term in self.terms:
        #    if term.prefactor != 0.0:
        #        acc.append(term)
        #        
        #self.terms = acc
            

    def gen_fortran_M1(self,si,sa,ref):
        code = '  ' + '! ### Spin case: i_'+si+', a_'+sa +' ###\n\n'
        for deltas,list_term in zip(self.unique_deltas,self.factorized_terms):
            d = []
            for delta in deltas:
                op1 = delta[0]
                op2 = delta[1]
                d.append([op1,op2])

            tmp = ''
            if d != []:
                tmp = str(d).replace('\'','')
                tmp = tmp.replace('[','(')
                tmp = tmp.replace(']',')')
                code  = '  !! Deltas:'+tmp+'\n'

            #code += '  !$OMP DO\n'
            shift = '  '
            code, shift = add_do_fortran(d,sa,'a',code,shift,False)
            code, shift = add_do_fortran(d,si,'i',code,shift,False)
            code = code + shift + self.which_M1(si,sa,d,ref)

            for term in list_term:
                p = str(term.prefactor)
                if p[0] != '-':
                    p = '+' + p
                p = p + 'd0'
                code = code + shift + p + term.Ts_to_fortran(shift) 

            code = code[:-4] + '\n'
            for i in range(len(shift)-2,0,-2):
                shft = ' '*i
                code = code + shft + 'enddo\n'
            #code += '  !$OMP ENDDO NOWAIT\n'

            print(code)

    def gen_fortran_M2_disconnected(self,si,sj,sa,sb,ref):   
        code = '  ' + '! ### Spin case: i_'+si+', j_'+sj+', a_'+sa+', b_'+sb+' ###\n\n'
        for deltas,list_term in zip(self.unique_deltas,self.factorized_terms):
            d = []
            for delta in deltas:
                op1 = delta[0]
                op2 = delta[1]
                d.append([op1,op2])
    
            tmp = ''
            if d != []:
                tmp = str(d).replace('\'','')
                tmp = tmp.replace('[','(')
                tmp = tmp.replace(']',')')
                code  = '  !! Deltas:'+tmp+'\n'
    
            #code += '  !$OMP DO\n'
            shift = '  '
            code, shift = add_do_fortran(d,sb,'b',code,shift,False)
            code, shift = add_do_fortran(d,sa,'a',code,shift,False)
            code, shift = add_do_fortran(d,sj,'j',code,shift,False)
            code, shift = add_do_fortran(d,si,'i',code,shift,False)
            code = code + shift + self.which_M2(si,sj,sa,sb,d,ref)
    
            for term in list_term:
                p = str(term.prefactor)
                if p[0] != '-':
                    p = '+' + p
                p = p + 'd0'
                disc = term.Ts[0] # I know that the disconnected term is in first position
                #print("disc:",disc.t)
                #if len(disc.t) != 2 or disc.t[0][0] == 'n' or disc.t[0][0] == 'm':
                #    print("That's not normal...)"+str(disc.t))
                #    sys.exit()
                #print('i'+si != disc.t[0] and 'j'+sj != disc.t[0],'i'+si ,'j'+sj, disc.t[0])
                #if 'i'+si != disc.t[0] and 'j'+sj != disc.t[0]:
                #    s1 = disc.t[0][1]
                #    if s1 == "a":
                #        l1 = "ma"
                #        l2 = "na"
                #    else:
                #        l1 = "nb"
                #        l2 = "mb"
                #    if 'a'+sa != disc.t[1] and 'b'+sb != disc.t[1]:
                #        str_term = shift + "if ("+str(disc.t[0])+" /= "+l1+" .or "+str(disc.t[1])+" /= "+l2+") then \n"
                #        str_term = str_term + "  " + shift + self.which_M2(si,sj,sa,sb,d,ref) + shift + "  " + p + term.Ts_to_fortran(shift)
                #        str_term = str_term[:-4] + " \n"
                #        str_term = str_term + shift + "endif   \n"
                #    else:
                #    #    print("Whaaaat???")
                #    #    sys.exit()
                #        str_term = shift + self.which_M2(si,sj,sa,sb,d,ref) + shift + p + term.Ts_to_fortran(shift)
                #        str_term = str_term[:-4] + "    \n"
                #else:
                #    str_term = shift + self.which_M2(si,sj,sa,sb,d,ref) + shift + p + term.Ts_to_fortran(shift)
                #    str_term = str_term[:-4] + "    \n"
                    
                code = code + shift + p + term.Ts_to_fortran(shift) 
    
            code = code[:-4] + '\n'
            for i in range(len(shift)-2,0,-2):
                shft = ' '*i
                code = code + shft + 'enddo\n'
            #code += '  !$OMP ENDDO NOWAIT\n'
    
            print(code)

    def gen_fortran_M2(self,si,sj,sa,sb,ref):   
        code = '  ' + '! ### Spin case: i_'+si+', j_'+sj+', a_'+sa+', b_'+sb+' ###\n\n'
        for deltas,list_term in zip(self.unique_deltas,self.factorized_terms):
            d = []
            for delta in deltas:
                op1 = delta[0]
                op2 = delta[1]
                d.append([op1,op2])
    
            tmp = ''
            if d != []:
                tmp = str(d).replace('\'','')
                tmp = tmp.replace('[','(')
                tmp = tmp.replace(']',')')
                code  = '  !! Deltas:'+tmp+'\n'
    
            #code += '  !$OMP DO\n'
            shift = '  '
            code, shift = add_do_fortran(d,sb,'b',code,shift,False)
            code, shift = add_do_fortran(d,sa,'a',code,shift,False)
            code, shift = add_do_fortran(d,sj,'j',code,shift,False)
            code, shift = add_do_fortran(d,si,'i',code,shift,False)
            code = code + shift + self.which_M2(si,sj,sa,sb,d,ref)
    
            for term in list_term:
                p = str(term.prefactor)
                if p[0] != '-':
                    p = '+' + p
                p = p + 'd0'
                code = code + shift + p + term.Ts_to_fortran(shift) 
    
            code = code[:-4] + '\n'
            for i in range(len(shift)-2,0,-2):
                shft = ' '*i
                code = code + shft + 'enddo\n'
            #code += '  !$OMP ENDDO NOWAIT\n'
    
            print(code)

    def which_M1(self,si,sa,d,ref):
        i = 'i'+si
        a = 'a'+sa
        M1 = 'M1_'+ref+'('+i+','+a+') = M1_'+ref+'('+i+','+a+') & \n'
        label = [i,a]
        for l in label:
            for elem in d:
                if l == elem[0]:
                    M1 = M1.replace(l,elem[1])
                elif l == elem[1]:
                    M1 = M1.replace(l,elem[0])
                
        return M1

    def which_M2(self,si,sj,sa,sb,d,ref):
        i = 'i'+si
        j = 'j'+sj
        a = 'a'+sa
        b = 'b'+sb
        M2 = 'M2_'+ref+'('+i+','+j+','+a+','+b+') = M2_'+ref+'('+i+','+j+','+a+','+b+') & \n'
        label = [i,j,a,b]
        for l in label:
            for elem in d:
                if l == elem[0]:
                    M2 = M2.replace(l,elem[1])
                elif l == elem[1]:
                    M2 = M2.replace(l,elem[0])
                
        return M2

def add_do_fortran(d,s_label,label,code,shift,disconnected):
    is_in = False
    for elem in d:
        #print(elem)
        for op in elem:
            #print(op,label+s_label)
            if label+s_label == op:
                is_in = True
    if not is_in:
        code = code + shift + 'do '+label+s_label+' = i_'+label+s_label+', f_'+label+s_label+'\n'
        shift = shift + '  '
                
        if label+s_label == 'ia' or label+s_label == 'ja':
            l = 'ma'
        elif label+s_label == 'aa' or label+s_label == 'ba':
            l = 'na'
        elif label+s_label == 'ib' or label+s_label == 'jb':
            l = 'nb'
        elif label+s_label == 'ab' or label+s_label == 'bb':
            l = 'mb'
        else:
            print('ooops')
            sys.exit()
            
        if s_label == 'a':
            if label == 'i' or label == 'j':
                txt = label+'b = '+label+'a + cc_nOa'
            elif label == 'a' or label == 'b':
                txt = label+'b = '+label+'a + cc_nVa'
            else:
                print('no such s_label')
                sys.exit()
                
        if s_label == 'b':
            if label == 'i' or label == 'j':
                txt = label+'a = '+label+'b - cc_nOa'
            elif label == 'a' or label == 'b':
                txt = label+'a = '+label+'b - cc_nVa'
            else:
                print('no such s_label')
                sys.exit()
        if s_label != 'a' and s_label != 'b':
            print('Well, we have a problem here')
            sys.exit()

        if not disconnected:
            code = code + shift + 'if ('+label+s_label+' == '+l+') cycle \n'
        code = code + shift + txt + '\n'
            
    return code, shift

def is_conflict_deltas(delta1,delta2):
    is_conflict = False
    for op1 in delta1:
        for op2 in delta2:
            #print(op1,op2,op1==op2)
            if op1 == op2:
                is_conflict = True
            
            #if op1[0] == op2[0] or op1[1] == op2[0]:
            #    count = count + 1
            #if op1[0] == op2[1] or op1[1] == op2[1]:
            #    count = count + 1
                
            #if count != 0:
            #    is_conflict = True
    #print(delta1,delta2,is_conflict)

    return is_conflict

def delta4_to_delta2(delta):
    d1 = delta[0][1] + delta[0][3]
    d2 = delta[1][1] + delta[1][3]
    d = [d1,d2]

    return d

def display_tex(tex):
    display(Latex(f'${tex}$'))

def apply_ops_eT(ops,Ts):
    if type(ops) != type(['aa','bb']):
        print('Error type arg ops in apply_ops_eT')
        sys.exit()
    if type(Ts) != type([['t1','t1','t2']]):
        print('Error type arg Ts in apply_ops_eT')
        sys.exit()
    
    # Init
    pq = pdaggerq.pq_helper("fermi")
    
    op_str = gen_left_str(ops)
    #print('Left ops:',op_str)
    
    # Set left operators
    pq.set_left_operators([[op_str]])
    #pq.set_left_operators([['e3(ira,isb,iig,aqa,apb,aag)']])
    
    #print('If there are many T, set the prefactor to 1/k! ...\n')
    # Set Ts operators
    #Ts = ['t1','t2']
    prefactor = 1.0/factorial(Ts.count('t1')) * 1.0/factorial(Ts.count('t2'))
    pq.add_operator_product(prefactor,Ts)
    #print(prefactor)
    #pq.add_operator_product(1.0/factorial(len(Ts)),Ts)
    #pq.add_operator_product(1.0,['t1','t1'])
    
    pq.simplify()
    
    # list of fully-contracted strings, then print
    terms = pq.fully_contracted_strings()
    #print(1,terms)
    #for term in terms:
    #    print(term)
    #    #pq.clear()
    #    obj = Term_pq(term)
    #    #print('prefactor',obj.prefactor)
    #    #print('T:',obj.l_T)
    #    obj.to_latex()

    return terms

def gen_all_T(max_rank,nb_min_op,nb_max_op):
    if type(max_rank) != type(1):
        print('Error type arg max_rank in gen_all_T')
        sys.exit()
    if type(nb_min_op) != type(1):
        print('Error type arg nb_min_op in gen_all_T')
        sys.exit()
    if type(nb_max_op) != type(1):
        print('Error type arg nb_max_op in gen_all_T')
        sys.exit()
        
    T = [[]]
    for i in range(1,max_rank+1):
        T[0].append([i])    
    #print(T)
    
    for j in range(1,nb_max_op//2+1):
        l = copy.deepcopy(T[j-1])
        idx = [i for i in range(1,max_rank+1)]
        res = []
        for elem in l:
            #print('e',elem)
            for i in idx:
                #print(sum(elem)+i)
                if i < elem[-1] or sum(elem)+i > nb_max_op//2: continue
                tmp = copy.deepcopy(elem)
                tmp.append(i)
                res.append(tmp)
            
        #print(res)
        if len(res) > 0: 
            T.append(res)
        #print(T)
        
    #print('T',T)
    #for elem in T:
    #    print(elem)
    
    # Reduce the number of dim
    acc = []
    for l in T:
        for ts in l:
            #print(ts)
            if sum(ts) < nb_min_op//2:
                continue
            acc.append(ts)
    T = acc

    # Transform to strings
    for i in range(len(T)):
        for j in range(len(T[i])):
            T[i][j] = 't'+str(T[i][j])

    return T

def spin_flip(list_op,idx_spin):
    if type(list_op) != type(['aa','bb']):
        print('Error type arg list_op in spin_flip')
        sys.exit()
    if type(idx_spin) != type(0):
        print('Error type arg idx_spin in spin_flip')
        sys.exit()
    
    acc = []
    for i in range(len(list_op)):
        tmp1 = copy.deepcopy(list_op[i][:idx_spin])
        s = list_op[i][idx_spin]
        if len(list_op[i]) > idx_spin:
            tmp2 = copy.deepcopy(list_op[i][idx_spin+1:])
        else:
            tmp2 = ''
            
        if s == 'b':
            s = 'a'
        elif s == 'a':
            s = 'b'
        acc.append(tmp1+s+tmp2)

    return acc

# Apply a spin flip on the result
def spin_flip_ltup(l,idx_perm_Ts):
    if type(l) != type([(1.0,[['aa']],[['aa']])]):
        print('Error type arg l in perm_ltup')
        sys.exit()
    if type(idx_perm_Ts) != type(1):
        print('Error type arg idx_perm_Ts in perm_ltup')
        sys.exit()

    list_prefactor = []
    list_deltas = []
    list_Ts = []
    l2 = copy.deepcopy(l)
    for prefactor, deltas, Ts in l2:
        list_prefactor.append(prefactor)
        list_deltas.append(deltas)
        list_Ts.append(Ts)

    acc_Ts = copy.deepcopy(list_Ts)
    for i in range(len(list_Ts)):
        for j in range(len(list_Ts[i])):
            #print(list_Ts[i][j])
            for k in range(len(list_Ts[i][j])):
                tmp1 = list_Ts[i][j][k][0:idx_perm_Ts]
                tmp2 = list_Ts[i][j][k][idx_perm_Ts+1:]
                tmp = list_Ts[i][j][k][idx_perm_Ts]
                if tmp == 'a':
                    tmp = 'b'
                elif tmp == 'b':
                    tmp = 'a'
                acc_Ts[i][j][k] = tmp1 + tmp + tmp2
            #print(list_Ts[i][j])

    acc_d = copy.deepcopy(list_deltas)
    for i in range(len(list_deltas)):
        for j in range(len(list_deltas[i])):
            for k in range(len(list_deltas[i][j])):
                tmp = list_deltas[i][j][k][0]
                tmp1 = list_deltas[i][j][k][1]
                tmp2 = list_deltas[i][j][k][3:]
                #print('1',list_deltas[i][j][k])
                if tmp == 'i':
                    tmp = 'a'
                    if k == 0:
                        op = '-'
                    else:
                        op = '+'
                elif tmp == 'a':
                    tmp = 'i'
                    if k == 0:
                        op = '+'
                    else:
                        op = '-'
                acc_d[i][j][k] = tmp + tmp1 + op + tmp2
                #print('2',list_deltas[i][j][k])
            #print('4',list_deltas[i][j])
                    
    acc = []
    for prefactor, deltas, Ts in zip(list_prefactor,acc_d,acc_Ts):
        acc.append((prefactor,deltas,Ts))

    return acc

# Search the unique elem of a list and sort them depending on their length
def search_unique_deltas(l):
    list_deltas = copy.deepcopy(l)
    if type(list_deltas) != type([['aa','bb'],['aa','bb']]):
        print('Error type arg list_deltas in search_unique_deltas')
        sys.exit()
    
    list_unique = []
    for deltas in list_deltas:
        if not is_in(deltas,list_unique):
            list_unique.append(deltas)

    # Sort
    ## Max len
    max_len = 0
    for elem in list_unique:
        if len(elem) > max_len:
            max_len = len(elem)

    ## Split depending on the length
    acc = [[] for i in range(max_len+1)]
    for deltas in list_unique:
        acc[len(deltas)].append(deltas)

    ## Reduction of the number of dimensions
    tmp = []
    for list_elem in acc:
        for elem in list_elem:
            tmp.append(elem)

    list_unique = tmp

    return list_unique

# Search if an elem is in a list l
def is_in(elem,l):
    for i in l:
        if elem == i:
            return True
    return False

def search_idx(elem,l):
    idx = 0
    for i in l:
        if i == elem:
            return idx
        idx = idx + 1

    print('Not found in the list')
    sys.exit()

def prod_ltup(l1,l2):
    if type(l1) != type([(1.0,[['aa','bb']],[['aa','bb']])]):
        print('Error type arg l1 in prod_ltup')
        sys.exit()
    if type(l2) != type(l1):
        print('Error type arg l1 in prod_ltup')
        sys.exit()
        
    # Prod
    res = []
    for prefactor1, deltas1, Ts1 in l1:
        for prefactor2, deltas2, Ts2 in l2:
            #print('1:',prefactor1,prefactor2)
            #print('1:',deltas1,deltas2)
            #print('1:',Ts1,Ts2)
            prefactor = prefactor1 * prefactor2
            deltas = copy.deepcopy(deltas1)
            is_in = False
            for elem in deltas2:
                for op2 in elem:
                    for e1 in deltas1:
                        for op1 in e1:
                            if op1 == op2:
                                is_in = True
                deltas.append(elem)
            # Remove the products of deltas leading to conflicts
            if is_in:
                continue
            deltas = remove_duplicate(deltas)
            Ts = copy.deepcopy(Ts1)
            for elem in Ts2:
                Ts.append(elem)
            Ts = sort_Ts(Ts)
            
            #print('2:',prefactor)
            #print('2:',deltas)
            #print('2:',Ts)
            res.append((prefactor,deltas,Ts))

    return res

def print_ltup(l):
    if type(l) != type([(1.0,[['aa','bb']],[['aa','bb']])]):
        print('Error type arg l in print_ltup')
        sys.exit()
        
    tex = ''
    for prefactor,deltas,Ts in l:
        tex = tex + prefactor_to_tex(prefactor) + ' \ '
        tex = tex + deltas_to_tex(deltas)
        tex = tex + Ts_to_tex(Ts)
    display(Latex(f'${tex}$'))

def print_fact_ltup(l,l_factor):
    if type(l) != type([(1.0,[['aa','bb']])]):
        print('Error type arg l in print_fact_ltup')
        sys.exit()
    if type(l_factor) != type([['aa','bb']]):
        print('Error type arg l_factor in print_fact_ltup')
        sys.exit()
        
    tex = '\\textbf{Factorized form:}'
    display(Latex(f'${tex}$'))
    for i in range(len(l)):
        tex = deltas_to_tex(l_factor[i]) + '\\bigl['
        #print('\nFactor:',list_unique_deltas[i])
        #if len(tex) > 0:
        #    display(Latex(f'${tex}$'))
        j = 0
        for elem in l[i]:
            #print(elem)
            e1,e2 = elem
            prefactor = prefactor_to_tex(e1)
            if j == 0 and prefactor[0] == '+':
                prefactor = prefactor[1:]
            tex = tex + prefactor + ' \ '
            tex = tex + Ts_to_tex(e2)
            j = j + 1
        tex = tex  + '\\bigr]'
        display(Latex(f'${tex}$'))

def print_contrib_M1(si,sa,len_res):
    if si == 'b':
        i = '\\bar{i}'
    else:
        i = 'i'
    if sa == 'b':
        a = '\\bar{a}'
    else:
        a = 'a'
        
    tex = 'M_{'+i+'}^{'+a+'} \\leftarrow '
    
    if len_res == 0:
        tex = tex + '0'
    null = display_tex(tex)
    #print(tex)
    
def print_contrib_M2(si,sj,sa,sb,len_res):
    if si == 'b':
        i = '\\bar{i}'
    else:
        i = 'i'
    if sj == 'b':
        j = '\\bar{j}'
    else:
        j = 'j'
    if sa == 'b':
        a = '\\bar{a}'
    else:
        a = 'a'
    if sb == 'b':
        b = '\\bar{b}'
    else:
        b = 'b'
    
    tex = 'M_{'+i+j+'}^{'+a+b+'} \\leftarrow '

    if len_res == 0:
        tex = tex + '0'
    null = display_tex(tex)
    #print(tex)

def factorize_from_ltup(l):
    if type(l) != type([(1.0,[['aa','bb']],[['aa','bb']])]):
        print('Error type arg l in factorize_from_ltup')
        sys.exit()
        
    list_prefactor = []
    list_deltas = []
    list_Ts = []
    for prefactor,deltas,Ts in l:
        list_prefactor.append(prefactor)
        list_deltas.append(deltas)
        list_Ts.append(Ts)

    return factorize_by_deltas(list_prefactor,list_deltas,list_Ts)

def factorize_by_deltas(list_prefactor,list_deltas,list_Ts):
    if type(list_prefactor) != type([1.0,1.0]):
        print('Error type arg list_prefactor in factorize_by_deltas')
        sys.exit()
    if type(list_deltas) != type([['aa','bb']]):
        print('Error type arg list_deltas in factorize_by_deltas')
        sys.exit()
    if type(list_Ts) != type([['aa','bb']]):
        print('Error type arg list_Ts in factorize_by_deltas')
        sys.exit()
        
    list_unique_deltas = search_unique_deltas(list_deltas)
    #for unique_deltas in list_unique_deltas:
    #    print(unique_deltas)
        
    # factorization
    fact_terms = [[] for i in range(len(list_unique_deltas))]
    for prefactor,deltas,Ts in zip(list_prefactor,list_deltas,list_Ts):
        idx = search_idx(deltas,list_unique_deltas)
        fact_terms[idx].append((prefactor,Ts))
        
    #for i in range(len(fact_terms)):
    #    print('\nFactor:',list_unique_deltas[i])
    #    for elem in fact_terms[i]:
    #        print(elem)

    return fact_terms, list_unique_deltas

def simplify_by_deltas(fact_terms):
    acc = []
    for terms in fact_terms:
        list_prefactor = []
        list_Ts = []
        for i in range(0,len(terms)):
            a_i = terms[i]
            prefactor_i = a_i[0]
            Ts_i = a_i[1]

            is_in = False
            for j in range(0,len(list_Ts)):
                prefactor_j = list_prefactor[j]
                Ts_j = list_Ts[j]
                if Ts_i == Ts_j:
                    list_prefactor[j] = prefactor_i + prefactor_j
                    is_in = True
                    break
            if is_in:
                continue
                
            #if prefactor_i == 0:
            #    continue
            list_prefactor.append(prefactor_i)
            list_Ts.append(Ts_i)
            
        tmp = []
        for p,t in zip(list_prefactor,list_Ts):
            #if p == 0.0:
            #    continue
            tmp.append((p,t))
        acc.append(tmp)

    return acc

def ordered_by_t1(fact_terms):

    f = []
    for terms in fact_terms:
        list_unique_inact_t1 = []
        for p,Ts in terms:
            for t in Ts:
                if len(t) != 2:
                    break
                
                if t[0][0] != 'i' and t[0][0] != 'j':
                    continue
                if t[1][0] != 'a' and t[1][0] != 'b':
                    continue
                is_in = False
                for elem in list_unique_inact_t1:
                    if t == elem:
                        is_in = True
                        break
                if is_in:
                    continue
            
                list_unique_inact_t1.append(t)

        acc = [[] for i in range(len(list_unique_inact_t1)+1)]
        for p,Ts in terms:
            for t in Ts:
                i = 0
                found = False
                for elem in list_unique_inact_t1:
                    if elem[0][0] == t[0][0] and elem[1][0] == t[1][0]:
                        found = True
                        break
                    i = i+1
                if found:
                    break
            
            if not found:
                i = -1
            
            acc[i].append((p,Ts))

        #print(acc)
        tmp = []
        for term in acc:
            for elem in term:
                tmp.append(elem)
    
        f.append(tmp)
        
    return f

def sort_Ts(Ts):
    if type(Ts) != type([['aa','bb']]):
        print('Error type arg Ts in sort Ts')
        sys.exit()
        
    max_len = 0
    for elem in Ts:
        if len(elem) > max_len:
            max_len = len(elem)

    acc = [[] for i in range(max_len+1)]

    for elem in Ts:
        #print(len(elem),len(acc))
        acc[len(elem)].append(elem)

    tmp = []
    for l in acc:
        for elem in l:
            tmp.append(elem)
            
    return tmp

def perm_op(op,idx_perm,label1,label2):
    if type(op) != type('aa'):
        print('Error type arg op in perm_op')
        sys.exit()
    if type(idx_perm) != type(0):
        print('Error type arg idx_perm in perm_op')
        sys.exit()
    if type(label1) != type('a'):
        print('Error type arg label1 in perm_op')
        sys.exit()
    if type(label2) != type('a'):
        print('Error type arg label1 in perm_op')
        sys.exit()
        
    #print('op:',op)
    tmp1 = op[:idx_perm]
    label = op[idx_perm]
    tmp2 = op[idx_perm+1:]
    #print(tmp1,label,tmp2)
    if label == label1:
        label = label2
    elif label == label2:
        label = label1
        
    return tmp1+label+tmp2

def perm_string(string,idx_perm,label1,label2):
    if type(string) != type(['aa']):
        print('Error type arg string in perm_string')
        sys.exit()
        
    acc = []
    for op in string:
        acc.append(perm_op(op,idx_perm,label1,label2))
        
    return acc

def perm_list(list_string,idx_perm,label1,label2):
    if type(list_string) != type([['aa']]):
        print('Error type arg list_string in perm_list')
        sys.exit()
        
    acc = []
    for string in list_string:
        #print('string',string)
        acc.append(perm_string(string,idx_perm,label1,label2))

    return acc

def perm_list_list(l,idx_perm,label1,label2):
    if type(l) != type([[['aa']]]):
        print('Error type arg l in perm_list_list')
        sys.exit()

    acc = []
    for elem in l:
        acc.append(perm_list(elem,idx_perm,label1,label2))

    return acc

def perm_ltup(l,sign,idx_perm_deltas,idx_perm_Ts,*label_tuples):
    if type(l) != type([(1.0,[['aa']],[['aa']])]):
        print('Error type arg l in perm_ltup')
        sys.exit()
    if type(sign) != type(1):
        print('Error type arg sign in perm_ltup')
        sys.exit()
    if type(label_tuples) != type(('a','b')) and type(label_tuples) != type([('a','b')]) :
        print('Error type arg label_tuples in perm_ltup')
        sys.exit()        

    list_prefactor = []
    list_deltas = []
    list_Ts = []
    l2 = copy.deepcopy(l)
    for prefactor, deltas, Ts in l2:
        list_prefactor.append(prefactor*sign)
        list_deltas.append(deltas)
        list_Ts.append(Ts)

    for label in label_tuples:
        list_deltas = perm_list_list(list_deltas,idx_perm_deltas,label[0],label[1])
        list_Ts = perm_list_list(list_Ts,idx_perm_Ts,label[0],label[1])

    acc = []
    for prefactor, deltas, Ts in zip(list_prefactor,list_deltas,list_Ts):
        acc.append((prefactor,deltas,Ts))

    return acc

def remove_duplicate(l):
    acc = []
    for i in l:
        is_in = False
        for j in acc:
            if i == j:
                is_in = True
        if not is_in:
            acc.append(i)

    return acc

# Remove disconnected terms by looking if there is at least one
# active index in each t
def remove_disconnected(l):
    if type(l) != type([(1.0,[[]],[[]])]):
        print('Type error arg l in function remove disconnected')
        sys.exit()

    l2 = copy.deepcopy(l)
    acc = []
    for prefactor, deltas, Ts in l2:
        for t in Ts:
            is_co = False
            for idx in t:
                if idx[0] == 'n' or idx[0] == 'm':
                    is_co = True
                    break
            if not is_co:
                break
        if not is_co:
            continue
        acc.append((prefactor,deltas,Ts))

    return acc

def do_calc2(factor,max_rank_t,sq_strings,ref):
    if type(factor) != type(1.0) and type(factor) != type(1):
        print('Erro type arg factor in do_calc')
        sys.exit()
    if type(max_rank_t) != type(1):
        print('Error type arg max_rank in do_calc')
        sys.exit()
    if type(sq_strings) != type([[['aa'],['aa']]]):
        print('Error type arg list_s in do_calc')
        sys.exit()

    list_act_idx = ['m','n']
        
    lterms = LTerms()
    for elem in sq_strings:
        #print(elem)
        list_sign, list_deltas, list_string = w.do_wick(elem)
        
        list_left_ops = []
        for sign,deltas,string in zip(list_sign,list_deltas,list_string):
            obj = w.Wicked_str(sign,deltas,string)
            #obj.eq_show()
            obj.crea_to_left()
            sign_ordered = obj.sign
            left_ops = obj.ops
            #obj.tex_show()
            #obj.eq_show()
        
            #print(left_ops)
            if len(left_ops) == 0:
                continue
            
            list_Ts = gen_all_T(max_rank_t,len(left_ops),len(left_ops))
            for Ts in list_Ts:
                #print('Left ops:',left_ops)
                #print('Ts',Ts)
                res = apply_ops_eT(left_ops,Ts)
                #if len(res) != 0: print('Res:',res)
            
                for elem in res:
                    #print(term)
                    obj2 = Term_pq(elem)
                    if obj2.prefactor == 0.0: 
                        continue
                    #print(obj2.prefactor,sign_ordered)
                    obj2.prefactor = obj2.prefactor * sign_ordered * factor
                    #print(obj2.prefactor)
                    obj2.to_latex()
                    #obj2.tex_show()
                    #print('')
                    #tex = tex + obj2.tex + deltas_to_tex(deltas)
                    #print(tex)

                    # Conversion to delta with 2 idx ops
                    d = []
                    for delta in deltas:
                        d.append(delta4_to_delta2(delta))

                    # Conversion to T obj
                    Ts = []
                    for elem in obj2.l_T:
                        t = T(elem,ref,list_act_idx)
                        Ts.append(t)
                        
                    term = Term(d,obj2.prefactor,Ts)
                    lterms.append_Term(term)
                    
    return lterms

def monoK_eTL_L(factor,mono_K,L_ref,max_rank_t):
    sq_strings = gen_all_orb_class_mono(mono_K)
    sq_strings = gen_all_spin_mono(sq_strings)
    
    spin = ['a','b']
    res = [[[] for i in spin] for a in spin]
    for i in range(len(spin)):
        for a in range(len(spin)):
            res[i][a] = do_calc2(factor,max_rank_t,sq_strings[i][a],L_ref)

    return res
                    

def gen_all_orb_class_mono(left_ops):
    orb = ['i','a']
    list_left_ops = []
    for i in orb:
        for a in orb:
            tmp = copy.deepcopy(left_ops)
            for k in range(len(tmp)):
                tmp[k] = tmp[k].replace('@(oi)',i)
                tmp[k] = tmp[k].replace('@(oa)',a)
            list_left_ops.append(tmp)
            
    return list_left_ops

def gen_all_spin_mono(list_left_ops):
    spin = ['a','b']
    list_left_ops_spin = [[0 for i in spin] for a in spin]
    for i in range(len(spin)):
        si = spin[i]
        for a in range(len(spin)):
            sa = spin[a]
            tmp = copy.deepcopy(list_left_ops)
            for k in range(len(tmp)):
                for l in range(len(tmp[k])):
                    tmp[k][l] = tmp[k][l].replace('@(si)',si)
                    tmp[k][l] = tmp[k][l].replace('@(sa)',sa)
            list_left_ops_spin[i][a] = tmp

    return list_left_ops_spin

def doubleK_eTL_L(factor,double_K,L_ref,max_rank_t):
    sq_strings = gen_all_orb_class_double(double_K)
    sq_strings = gen_all_spin_double(sq_strings)

    spin = ['a','b']
    res = [[[[0 for i in spin] for j in spin] for a in spin] for b in spin]
    for i in range(len(spin)):
        for j in range(len(spin)):
            for a in range(len(spin)):
                for b in range(len(spin)):
                    res[i][j][a][b] = do_calc2(factor,max_rank_t,sq_strings[i][j][a][b],L_ref)

    return res
                    

def gen_all_orb_class_double(left_ops):
    orb = ['i','a']
    list_left_ops = []
    for i in orb:
        for j in orb:
            for a in orb:
                for b in orb:
                    tmp = copy.deepcopy(left_ops)
                    for k in range(len(tmp)):
                        tmp[k] = tmp[k].replace('@(oi)',i)
                        tmp[k] = tmp[k].replace('@(oj)',j)
                        tmp[k] = tmp[k].replace('@(oa)',a)
                        tmp[k] = tmp[k].replace('@(ob)',b)
                    list_left_ops.append(tmp)
    
    return list_left_ops

def gen_all_spin_double(list_left_ops):
    spin = ['a','b']
    list_left_ops_spin = [[[[0 for i in spin] for j in spin] for a in spin] for b in spin]
    for i in range(len(spin)):
        si = spin[i]
        for j in range(len(spin)):
            sj= spin[j]
            for a in range(len(spin)):
                sa = spin[a]
                for b in range(len(spin)):
                    sb = spin[b]
                    tmp = copy.deepcopy(list_left_ops)
                    for k in range(len(tmp)):
                        for l in range(len(tmp[k])):
                            tmp[k][l] = tmp[k][l].replace('@(si)',si)
                            tmp[k][l] = tmp[k][l].replace('@(sj)',sj)
                            tmp[k][l] = tmp[k][l].replace('@(sa)',sa)
                            tmp[k][l] = tmp[k][l].replace('@(sb)',sb)
                    list_left_ops_spin[i][j][a][b] = tmp
    return list_left_ops_spin
