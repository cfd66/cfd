#include "udf.h"

#define A_H2O 1e-5  /* 水蒸气渗透率 [mol/(m²·s·Pa)] */
#define H2O_index 0
#define UDM_FLUX_RCT 0  /* UDM存储反应侧通量值 */
#define UDM_FLUX_PERM 1 /* UDM存储渗透侧通量值 */

/*--------------------------------------------------------------*/
/*                    反应侧源项UDF (质量移除)                   */
/*--------------------------------------------------------------*/
DEFINE_SOURCE(H2O_source_reactant, c, t, dS, eqn)
{
    real src;
    real vol = C_VOLUME(c, t);
    
    /* 从UDM获取总质量损失率 [kg/s] */
    real mass_loss_rate = C_UDMI(c, t, UDM_FLUX_RCT);
    
    /* 计算源项 (负值，单位体积的质量损失率) [kg/(m³·s)] */
    src = -mass_loss_rate / vol;
    
    /* 导数项设为0 (源项不直接依赖于求解变量) */
    dS[eqn] = 0.0;  
    
    return src;
}

/*--------------------------------------------------------------*/
/*                    渗透侧源项UDF (质量增加)                   */
/*--------------------------------------------------------------*/
DEFINE_SOURCE(H2O_source_permeate, c, t, dS, eqn)
{
    real src;
    real vol = C_VOLUME(c, t);
    
    /* 从UDM获取总质量增加率 [kg/s] */
    real mass_gain_rate = C_UDMI(c, t, UDM_FLUX_PERM);
    
    /* 计算源项 (正值，单位体积的质量增加率) [kg/(m³·s)] */
    src = mass_gain_rate / vol;
    
    /* 导数项设为0 */
    dS[eqn] = 0.0;  
    
    return src;
}

/*--------------------------------------------------------------*/
/*              UDM初始化宏 (在每个迭代步开始时清零)            */
/*--------------------------------------------------------------*/
DEFINE_EXECUTE_AT_END(reset_flux_udms)
{
    Domain *domain=Get_Domain(1);
    Thread *t;
    cell_t c;
    
    /* 遍历所有计算域线程 */
    thread_loop_c(t, domain)
    {
        if(THREAD_STORAGE(t,SV_UDM_I)!=NULL){
            begin_c_loop(c, t)  // 修正缩进
            {
                C_UDMI(c, t, UDM_FLUX_RCT) = 0.0;
                C_UDMI(c, t, UDM_FLUX_PERM) = 0.0;
            }
            end_c_loop(c, t)  // 修正缩进
        }
    }
}

/*--------------------------------------------------------------*/
/*           修改后的边界通量UDF (添加UDM存储功能)               */
/*--------------------------------------------------------------*/
DEFINE_PROFILE(H2O_selective_membrane, thread, index)
{
    face_t f;
    cell_t c0, c1;
    Thread *thread0, *thread1;
    real p_H2O_ret, p_H2O_perm;
    real J_mol;
    
    thread0 = THREAD_T0(thread);  // 反应侧（CO2+H2O）
    thread1 = THREAD_T1(thread);  // 渗透侧（N2+H2O）
    
    real p_op = RP_Get_Real("operating-pressure");
    
    begin_f_loop(f, thread)
    {   
        c0 = F_C0(f, thread);
        c1 = F_C1(f, thread);
        
        if (!c0 || !c1) continue;

        int i = 0;
        real sum_mole_ret = 0.0;
        real sum_mole_perm = 0.0;
        Material *sp = NULL;
        real Mw_H2O = 0.0;
        real area[ND_ND];
        
        /* 反应侧分压计算 */
        mixture_species_loop(THREAD_MATERIAL(thread0), sp, i)
        {
            real Mwi = MATERIAL_PROP(sp, PROP_mwi);
            if(i == H2O_index) Mw_H2O = Mwi;
            sum_mole_ret += C_YI(c0, thread0, i) / Mwi;  
        }
        p_H2O_ret = (C_P(c0, thread0) + p_op) * 
                   (C_YI(c0, thread0, H2O_index) / Mw_H2O) / sum_mole_ret;
        
        /* 渗透侧分压计算 */
        mixture_species_loop(THREAD_MATERIAL(thread1), sp, i)
        {
            real Mwi = MATERIAL_PROP(sp, PROP_mwi);
            if(i == H2O_index) Mw_H2O = Mwi;
            sum_mole_perm += C_YI(c1, thread1, i) / Mwi;
        }
        p_H2O_perm = (C_P(c1, thread1) + p_op) * 
                    (C_YI(c1, thread1, H2O_index) / Mw_H2O) / sum_mole_perm;                 
        
        /* 单向渗透保护 */
        real delta_p = p_H2O_ret - p_H2O_perm;
        J_mol = (delta_p > 0) ? A_H2O * delta_p : 0.0;
        
        /* 将摩尔通量转换为质量通量 [kg/(m²·s)] */
        real J_kg = J_mol * (Mw_H2O / 1000.0);
        F_PROFILE(f, thread, index) = J_kg;
        
        /* 计算当前面的质量流量 [kg/s] */
        F_AREA(area, f, thread);  // 先调用F_AREA宏填充数组
        real face_area = NV_MAG(area);  // 单独计算向量模长
        real mass_flow = J_kg * face_area;
        
        /* 将质量流量分配到相邻单元的UDM中 */
        C_UDMI(c0, thread0, UDM_FLUX_RCT) += mass_flow;  // 反应侧损失质量
        C_UDMI(c1, thread1, UDM_FLUX_PERM) += mass_flow;  // 渗透侧获得质量
    }
    end_f_loop(f, thread)
}