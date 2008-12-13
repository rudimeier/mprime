// WorkerDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "WorkerDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// In theory, the maximum number of workers should be number of logical cpus.
// However, local.ini could specify more worker threads than logical cpus (for
// example, when local.ini is copied from a dual-core to a single-core machine).
// We must let the user manipulate the options on these worker threads that
// don't have a CPU to run on.

unsigned int max_num_workers (void)
{
	return (max (NUM_WORKER_THREADS, NUM_CPUS * CPU_HYPERTHREADS));
}

/////////////////////////////////////////////////////////////////////////////
// CWorkerDlg dialog


CWorkerDlg::CWorkerDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CWorkerDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CWorkerDlg)
	m_num_thread = 1;
	m_priority = 1;
	thread_num = 0;
	//}}AFX_DATA_INIT
}

void CWorkerDlg::DoDataExchange(CDataExchange* pDX)
{

	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CWorkerDlg)
	DDX_Control(pDX, IDC_THREAD_TEXT, c_num_thread_text);
	DDX_Control(pDX, IDC_THREAD, c_num_thread);
	DDX_Text(pDX, IDC_THREAD, m_num_thread);
	DDV_MinMaxUInt(pDX, m_num_thread, 1, max_num_workers ());
	DDX_Text(pDX, IDC_PRIORITY, m_priority);
	DDV_MinMaxUInt(pDX, m_priority, 1, 10);
	DDX_Control(pDX, IDC_WORKGROUP, c_workgroup);

	DDX_Control(pDX, IDC_THREADNUM_TEXT, c_threadnum_text);
	DDX_Control(pDX, IDC_COMBO1, c_threadnum);

	DDX_Control(pDX, IDC_WORK_TYPE_TEXT, c_work_pref_text);
	DDX_Control(pDX, IDC_COMBO2, c_work_pref);

	DDX_Control(pDX, IDC_AFFINITY_TEXT, c_affinity_text);
	DDX_Control(pDX, IDC_COMBO3, c_affinity);

	DDX_Control(pDX, IDC_NUMCPUS_TEXT, c_numcpus_text);
	DDX_Control(pDX, IDC_NUMCPUS, c_numcpus);

	DDX_Control(pDX, IDC_WARN1, c_warn1_text);
	DDX_Control(pDX, IDC_WARN2, c_warn2_text);
	DDX_Control(pDX, IDC_WARN3, c_warn3_text);

	//}}AFX_DATA_MAP
	c_num_thread_text.EnableWindow (max_num_workers () > 1);
	c_num_thread.EnableWindow (max_num_workers () > 1);
	c_workgroup.EnableWindow (USE_PRIMENET || NUM_CPUS * CPU_HYPERTHREADS > 1);
	c_threadnum_text.EnableWindow (m_num_thread > 1);
	c_threadnum.EnableWindow (m_num_thread > 1);
	c_work_pref_text.EnableWindow (USE_PRIMENET);
	c_work_pref.EnableWindow (USE_PRIMENET);
	c_affinity_text.EnableWindow (NUM_CPUS * CPU_HYPERTHREADS > 1);
	c_affinity.EnableWindow (NUM_CPUS * CPU_HYPERTHREADS > 1);
	c_numcpus_text.EnableWindow (m_num_thread < NUM_CPUS * CPU_HYPERTHREADS || !AreAllTheSame (m_numcpus) || m_numcpus[0] != 1);
	c_numcpus.EnableWindow (m_num_thread < NUM_CPUS * CPU_HYPERTHREADS || !AreAllTheSame (m_numcpus) || m_numcpus[0] != 1);
}


BEGIN_MESSAGE_MAP(CWorkerDlg, CDialog)
	//{{AFX_MSG_MAP(CWorkerDlg)
	ON_EN_KILLFOCUS(IDC_THREAD, &CWorkerDlg::OnEnKillfocusNumThread)
	ON_EN_SETFOCUS(IDC_PRIORITY, &CWorkerDlg::OnEnSetfocusPriority)
	ON_EN_KILLFOCUS(IDC_PRIORITY, &CWorkerDlg::OnEnKillfocusPriority)
	ON_CBN_SELCHANGE(IDC_COMBO1, &CWorkerDlg::OnCbnKillfocusThreadNum)
	ON_CBN_KILLFOCUS(IDC_COMBO1, &CWorkerDlg::OnCbnKillfocusThreadNum)
	ON_CBN_KILLFOCUS(IDC_COMBO2, &CWorkerDlg::OnCbnKillfocusWorkType)
	ON_CBN_KILLFOCUS(IDC_COMBO3, &CWorkerDlg::OnCbnKillfocusAffinity)
	ON_EN_KILLFOCUS(IDC_NUMCPUS, &CWorkerDlg::OnEnKillfocusNumCpus)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CWorkerDlg internal routines

int CWorkerDlg::AreAllTheSame (
	int	*array)
{
	int	i;

	for (i = 1; i < (int) m_num_thread; i++)
		if (array[i-1] != array[i]) return (FALSE);
	return (TRUE);
}

int map_work_pref_to_sel (
	int	work_pref)
{
	switch (work_pref) {
	case PRIMENET_WP_WHATEVER:
		return (0);
	case PRIMENET_WP_WORLD_RECORD:
		return (1);
	case PRIMENET_WP_LL_FIRST:
		return (2);
	case PRIMENET_WP_LL_DBLCHK:
		return (3);
	case PRIMENET_WP_FACTOR:
		return (4);
	case PRIMENET_WP_PFACTOR:
		return (5);
	case PRIMENET_WP_FACTOR_LMH:
		return (6);
	case PRIMENET_WP_ECM_SMALL:
		return (7);
	case PRIMENET_WP_ECM_FERMAT:
		return (8);
	case PRIMENET_WP_LL_100M:
		return (9);
	default:
		return (10);
	}
}

int map_sel_to_work_pref (
	int	sel)
{
	switch (sel) {
	case 0:
		return (PRIMENET_WP_WHATEVER);
	case 1:
		return (PRIMENET_WP_WORLD_RECORD);
	case 2:
		return (PRIMENET_WP_LL_FIRST);
	case 3:
		return (PRIMENET_WP_LL_DBLCHK);
	case 4:
		return (PRIMENET_WP_FACTOR);
	case 5:
		return (PRIMENET_WP_PFACTOR);
	case 6:
		return (PRIMENET_WP_FACTOR_LMH);
	case 7:
		return (PRIMENET_WP_ECM_SMALL);
	case 8:
		return (PRIMENET_WP_ECM_FERMAT);
	case 9:
		return (PRIMENET_WP_LL_100M);
	}
	return (-1);
}

void CWorkerDlg::InitComboBoxText (void)
{
	int	i, sel;
	char	buf[80];

// Populate the thread number combo box

	c_threadnum_text.EnableWindow (m_num_thread > 1);
	c_threadnum.EnableWindow (m_num_thread > 1);
	c_threadnum.ResetContent ();
	if (m_num_thread == 1) {
		thread_num = 1;
		c_threadnum.AddString ("Worker #1");
		c_threadnum.SetCurSel (0);
	} else {
		if (thread_num > (int) m_num_thread) thread_num = 0;
		c_threadnum.AddString ("All workers");
		for (i = 1; i <= (int) m_num_thread; i++) {
			sprintf (buf, "Worker #%d", i);
			c_threadnum.AddString (buf);
		}
		c_threadnum.SetCurSel (thread_num);
	}

// Populate the work type combo box

	c_work_pref.ResetContent ();
	if (thread_num)
		sel = map_work_pref_to_sel (m_work_pref[thread_num-1]);
	else if (AreAllTheSame (m_work_pref))
		sel = map_work_pref_to_sel (m_work_pref[0]);
	else {
		c_work_pref.AddString ("Mixed Settings");
		sel = 0;
	}
	c_work_pref.AddString ("Whatever makes the most sense");
	c_work_pref.AddString ("World record sized numbers to test");
	c_work_pref.AddString ("First time tests");
	c_work_pref.AddString ("Double-check tests");
	c_work_pref.AddString ("Trial factoring");
	c_work_pref.AddString ("P-1 factoring");
	c_work_pref.AddString ("Trial factoring to low limits");
	c_work_pref.AddString ("ECM on small Mersenne numbers");
	c_work_pref.AddString ("ECM on Fermat numbers");
	c_work_pref.AddString ("100,000,000 digit numbers to test");
	if (sel == 10) c_work_pref.AddString ("Other");
	c_work_pref.SetCurSel (sel);

// Populate the affinity combo box.  If there is only one worker thread,
// then smart affinity assignment is the same as run on any cpu.

	c_affinity.ResetContent ();
	if (thread_num) {
		sel = m_affinity[thread_num-1];
		if (m_num_thread > 1)
			sel = (sel == 99) ? 0 : (sel == 100) ? 1 : sel + 2;
		else
			sel = (sel == 99 || sel == 100) ? 0 : sel + 1;
	} else if (AreAllTheSame (m_affinity)) {
		sel = m_affinity[0];
		sel = (sel == 99) ? 0 : (sel == 100) ? 1 : sel + 2;
	} else {
		c_affinity.AddString ("Mixed Settings");
		sel = 0;
	}
	c_affinity.AddString ("Run on any CPU");
	if (m_num_thread > 1) c_affinity.AddString ("Smart assignment");
	if (thread_num) {
		for (i = 1; i <= (int) (NUM_CPUS * CPU_HYPERTHREADS); i++) {
			sprintf (buf, "CPU #%d", i);
			c_affinity.AddString (buf);
		}
	}
	c_affinity.SetCurSel (sel);

// Populate the num cpus edit box.

	if (thread_num) {
		sprintf (buf, "%d", m_numcpus[thread_num-1]);
		c_numcpus.SetWindowText (buf);
	} else if (AreAllTheSame (m_numcpus)) {
		sprintf (buf, "%d", m_numcpus[0]);
		c_numcpus.SetWindowText (buf);
	} else {
		c_numcpus.SetWindowText ("Mixed");
	}
	c_numcpus_text.EnableWindow (m_num_thread < NUM_CPUS * CPU_HYPERTHREADS || !AreAllTheSame (m_numcpus) || m_numcpus[0] != 1);
	c_numcpus.EnableWindow (m_num_thread < NUM_CPUS * CPU_HYPERTHREADS || !AreAllTheSame (m_numcpus) || m_numcpus[0] != 1);
}


/////////////////////////////////////////////////////////////////////////////
// CWorkerDlg message handlers

BOOL CWorkerDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	c_warn1_text.EnableWindow (FALSE);
	c_warn2_text.EnableWindow (FALSE);
	c_warn3_text.EnableWindow (FALSE);

	InitComboBoxText ();

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CWorkerDlg::OnEnKillfocusNumThread()
{
	char	buf[80];
	unsigned int num_thread;

	c_num_thread.GetWindowText (buf, sizeof (buf));
	num_thread = atoi (buf);
	if (num_thread < 1) num_thread = 1;
	if (num_thread > max_num_workers ())
		num_thread = max_num_workers ();
	sprintf (buf, "%d", num_thread);
	c_num_thread.SetWindowText (buf);
	if (num_thread != m_num_thread) {
		m_num_thread = num_thread;
		InitComboBoxText ();
				// Might change thread_num if we reduced
				// number of worker threads
	}
}

void CWorkerDlg::OnEnSetfocusPriority()
{
	c_warn1_text.EnableWindow (TRUE);
	c_warn2_text.EnableWindow (TRUE);
	c_warn3_text.EnableWindow (TRUE);
}

void CWorkerDlg::OnEnKillfocusPriority()
{
	c_warn1_text.EnableWindow (FALSE);
	c_warn2_text.EnableWindow (FALSE);
	c_warn3_text.EnableWindow (FALSE);
}

void CWorkerDlg::OnCbnKillfocusThreadNum()
{
	thread_num = c_threadnum.GetCurSel ();
	InitComboBoxText ();	// Display data for entered thread num
}

void CWorkerDlg::OnCbnKillfocusWorkType()
{
	int	sel, i, work_pref;

	sel = c_work_pref.GetCurSel ();
	if (thread_num) {
		work_pref = map_sel_to_work_pref (sel);
		if (work_pref != -1) m_work_pref[thread_num-1] = work_pref;
	}
	else if (AreAllTheSame (m_work_pref)) {
		work_pref = map_sel_to_work_pref (sel);
		if (work_pref != -1)
			for (i = 0; i < MAX_NUM_WORKER_THREADS; i++)
				m_work_pref[i] = work_pref;
	} else if (sel) {
		work_pref = map_sel_to_work_pref (sel-1);
		for (i = 0; i < MAX_NUM_WORKER_THREADS; i++)
			m_work_pref[i] = work_pref;
	}
}

void CWorkerDlg::OnCbnKillfocusAffinity()
{
	int	sel, i, affinity;

	sel = c_affinity.GetCurSel ();
	if (thread_num) {
		if (m_num_thread > 1)
			affinity = (sel == 0) ? 99 : (sel == 1) ? 100 : sel - 2;
		else
			affinity = (sel >= 1) ? sel - 1 :
				   (m_affinity[0] == 100) ? 100 : 99;
		m_affinity[thread_num-1] = affinity;
	} else if (AreAllTheSame (m_affinity)) {
		affinity = (sel == 0) ? 99 : (sel == 1) ? 100 : sel - 2;
		for (i = 0; i < MAX_NUM_WORKER_THREADS; i++)
			m_affinity[i] = affinity;
	} else if (sel) {
		affinity = (sel == 1) ? 99 : (sel == 2) ? 100 : sel - 3;
		for (i = 0; i < MAX_NUM_WORKER_THREADS; i++)
			m_affinity[i] = affinity;
	}
}

void CWorkerDlg::OnEnKillfocusNumCpus()
{
	char	buf[80];
	unsigned int num_cpus;
	int	i;

	c_numcpus.GetWindowText (buf, sizeof (buf));
	if (buf[0] >= '0' && buf[0] <= '9') {
		num_cpus = atoi (buf);
		if (num_cpus < 1) num_cpus = 1;
		if (num_cpus > NUM_CPUS * CPU_HYPERTHREADS)
			num_cpus = NUM_CPUS * CPU_HYPERTHREADS;
		sprintf (buf, "%d", num_cpus);
		c_numcpus.SetWindowText (buf);
		if (thread_num)
			m_numcpus[thread_num-1] = num_cpus;
		else
			for (i = 0; i < MAX_NUM_WORKER_THREADS; i++)
				m_numcpus[i] = num_cpus;
	}
	c_numcpus_text.EnableWindow (m_num_thread < NUM_CPUS * CPU_HYPERTHREADS || !AreAllTheSame (m_numcpus) || m_numcpus[0] != 1);
	c_numcpus.EnableWindow (m_num_thread < NUM_CPUS * CPU_HYPERTHREADS || !AreAllTheSame (m_numcpus) || m_numcpus[0] != 1);
}

