import numpy as _np
import pykep as pk

def primer_vector(DVi, DVj, Mji, Mjk):
    """This function computes the primer vector in a point k, relative
    to finite impulses in i and j.
    
    Args:
        *DVi* (:class:`ndarray` - (3,)): the impulse at point i.
        
        *DVj* (:class:`ndarray` - (3,)): the impulse at point j.
        
        *Mji* (:class:`ndarray` - (6,6)): the state transition matrix from i to j (dxj = Mji dxi).
        
        *Mjk* (:class:`ndarray` - (6,6)): the state transition matrix from k to j (dxj = Mji dxi).
        
        Returns:
            :class:`tuple`: The primer vector, the Aik matrix, the Ajk matrix.
            
        Note:
            The impulse transfer matrix Anm is defined as those matrices that allow to compute the
            variation of the impulse at point n given the variation of the impulse at point m. In formal terms, 
            dDVn = Anm dDVm. All variations are such that the terminal state is kept fixed.
    """
    Aik = -(_np.linalg.inv(Mji[:3,3:]))@Mjk[:3,3:]
    Ajk = -(Mji[3:,3:]@Aik + Mjk[3:,3:])
    p = - Aik.T@DVi/_np.linalg.norm(DVi) - Ajk.T@DVj/_np.linalg.norm(DVj)
    return p, Aik, Ajk

def primer_vector_surrogate(DVk, Mki, Mkj):
    """This function computes the surrogate primer vector at two points i and j, 
    corresponding to a single finite impulse at point k.
    
    Args:
        *DVk* (:class:`ndarray` - (3,)): the impulse at point k.
        
        *Mki* (:class:`ndarray` - (6,6)): the state transition matrix from i to k (dxk = Mki dxi).
        
        *Mkj* (:class:`ndarray` - (6,6)): the state transition matrix from j to k (dxk = Mkj dxj).
        
        Returns:
            :class:`tuple`: The surrogate primer vector, the Aij matrix, the Akj matrix.
            
        Note:
            The impulse transfer matrix Anm is defined as that matrix that allows to compute 
            (in this surrogate case) the variation of the impulse at point n given the 
            variation of the impulse at point m. In formal terms, dDVn = Anm dDVm. 
            All variations are such that the terminal state is kept fixed.
    """
    Aij = -(_np.linalg.inv(Mki[:3, 3:]) @ Mkj[:3, 3:])
    Akj = -(Mki[3:, 3:] @ Aij + Mkj[3:, 3:])
    B = Aij
    b =  - Akj.T@DVk / _np.linalg.norm(DVk)
    p_surrogate_norm, u_star = pk.trajopt.minBu_bu_p(B, b)
    p_surrogate = p_surrogate_norm * u_star
    return p_surrogate, Aij, Akj
    
