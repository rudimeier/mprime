#pragma once
#include "afxwin.h"


// CTortureDlg dialog

class CTortureDlg : public CDialog
{
	DECLARE_DYNAMIC(CTortureDlg)

public:
	CTortureDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~CTortureDlg();

// Dialog Data
	enum { IDD = IDD_TORTURE };
	CStatic	c_thread_text;
	CEdit	c_thread;
	UINT	m_thread;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	int m_torture_type;
	int m_minfft;
	int m_maxfft;
	BOOL m_in_place_fft;
	int m_memory;
	int m_timefft;
	CStatic c_minfft_text;
	CEdit c_minfft;
	CStatic c_maxfft_text;
	CEdit c_maxfft;
	CButton c_in_place_fft;
	CStatic c_memory_text;
	CEdit c_memory;
	CStatic c_timefft_text;
	CEdit c_timefft;

	int	m_blendmemory;

	afx_msg void OnBnClickedL2Cache();
	afx_msg void OnBnClickedInPlace();
	afx_msg void OnBnClickedBlend();
	afx_msg void OnBnClickedCustom();
	afx_msg void OnBnClickedInPlaceFft();
};
